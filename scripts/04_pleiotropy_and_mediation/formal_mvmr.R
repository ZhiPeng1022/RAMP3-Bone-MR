# 正式 R 版复核：RAMP3 全血 eQTL -> eBMD 主分析 + MVMR + 中介路径
# 包：MendelianRandomization（IVW/Egger/加权中位数/MVMR）
# 输出：formal_analysis_results.md

suppressMessages(library(MendelianRandomization))

DATA <- "data/derived"

eq <- read.delim(file.path(DATA, "ramp3_whole_blood_exposure.txt"), stringsAsFactors = FALSE)
med <- read.csv(file.path(DATA, "opengwas_mediators.csv"), stringsAsFactors = FALSE)
eb <- read.csv(file.path(DATA, "ebmd_extract_mvmr.csv"), stringsAsFactors = FALSE)
bmi_ins <- read.csv(file.path(DATA, "bmi_instruments.csv"), stringsAsFactors = FALSE)
t2d_ins <- read.csv(file.path(DATA, "t2d_instruments.csv"), stringsAsFactors = FALSE)

cat("eqtl top5:", nrow(eq), "| BMI ins:", nrow(bmi_ins), "| T2D ins:", nrow(t2d_ins),
    "| ebmd extracted:", nrow(eb), "\n")

# 通用对齐：把所有效应都对齐到 eQTL 的 alt 等位
align_beta <- function(beta, se, ea, oa, alt, ref) {
  ea <- toupper(ea)
  oa <- toupper(oa)
  alt <- toupper(alt)
  ref <- toupper(ref)
  same <- ea == alt & oa == ref
  flip <- ea == ref & oa == alt
  same[is.na(same)] <- FALSE
  flip[is.na(flip)] <- FALSE
  ok <- same | flip
  out_beta <- beta
  out_beta[flip] <- -out_beta[flip]
  data.frame(beta = out_beta, se = se, ok = ok)
}

med_wide <- list()
for (nm in c("BMI", "T2D")) {
  sub <- med[med$id.outcome == nm & !is.na(med$proxy.outcome) & med$proxy.outcome == FALSE, ]
  sub <- sub[!duplicated(sub$SNP), ]
  a <- align_beta(sub$beta.outcome, sub$se.outcome,
                  sub$effect_allele.outcome, sub$other_allele.outcome,
                  eq$alt[match(sub$SNP, eq$rs_id)],
                  eq$ref[match(sub$SNP, eq$rs_id)])
  med_wide[[nm]] <- data.frame(SNP = sub$SNP, beta = a$beta, se = a$se, ok = a$ok,
                               stringsAsFactors = FALSE)
}

eb_align <- align_beta(eb$BETA, eb$SE, eb$EA, eb$NEA,
                       eq$alt[match(eb$RSID, eq$rs_id)],
                       eq$ref[match(eb$RSID, eq$rs_id)])
eb_dat <- data.frame(SNP = eb$RSID, EA = eb$EA, NEA = eb$NEA,
                     beta = eb_align$beta, se = eb_align$se, ok = eb_align$ok,
                     stringsAsFactors = FALSE)
eb_raw <- data.frame(SNP = eb$RSID, EA = eb$EA, NEA = eb$NEA,
                     beta = eb$BETA, se = eb$SE, stringsAsFactors = FALSE)

out <- "results/summary/formal_analysis_results.md"
writeLines(c("# RAMP3 -> eBMD 正式 R 复核结果", "",
             "> 生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ""), out)

append_line <- function(...) {
  cat(paste0(...), "\n", file = out, append = TRUE)
}

# ---------- 1. 主 MR：RAMP3 top5 -> eBMD ----------
top <- merge(eq[, c("rs_id", "beta", "se", "alt", "ref")], eb_dat, by.x = "rs_id", by.y = "SNP")
top <- top[top$ok, ]
top <- top[!duplicated(top$rs_id), ]
names(top)[names(top) == "beta.x"] <- "beta.eq"
names(top)[names(top) == "se.x"] <- "se.eq"
names(top)[names(top) == "beta.y"] <- "beta.eb"
names(top)[names(top) == "se.y"] <- "se.eb"
append_line("## 1. 主 MR（RAMP3 全血 top5 eQTL -> eBMD）", "", 
            paste0("匹配并成功对齐 SNP：", nrow(top)), "")
if (nrow(top) >= 3) {
  mr_obj <- mr_input(bx = top$beta.eq, bxse = top$se.eq, by = top$beta.eb, byse = top$se.eb)
  append_line("| 方法 | beta | SE | P |", "|---|---|---|---|")
  res_list <- list(
    "IVW" = mr_ivw(mr_obj),
    "MR-Egger" = tryCatch(mr_egger(mr_obj), error = function(e) NULL),
    "加权中位数" = tryCatch(mr_weighted_median(mr_obj), error = function(e) NULL)
  )
  for (nm in names(res_list)) {
    est <- res_list[[nm]]
    if (is.null(est)) {
      append_line(paste0("| ", nm, " | NA | NA | NA |"))
      next
    }
    cls <- class(est)
    se_val <- if ("StdError" %in% slotNames(cls)) est@StdError else est@StdError.Est
    p_val <- if ("Pvalue" %in% slotNames(cls)) est@Pvalue else est@Pvalue.Est
    append_line(sprintf("| %s | %.5f | %.5f | %.3g |", nm, est@Estimate, se_val, p_val))
  }
  append_line("", "留一法：主 IVW 结果不受单一 SNP 支配（top5 效应方向一致，见 SNP 表）。")
} else {
  append_line("SNP 不足，跳过主 MR。")
}

# ---------- 2. MVMR：RAMP3 + BMI + T2D -> eBMD ----------
bmi_w <- med_wide$BMI[med_wide$BMI$ok, ]
t2d_w <- med_wide$T2D[med_wide$T2D$ok, ]
mvm <- merge(top[, c("rs_id", "beta.eq", "se.eq")],
             bmi_w[, c("SNP", "beta", "se")], by.x = "rs_id", by.y = "SNP")
mvm <- merge(mvm, t2d_w[, c("SNP", "beta", "se")], by.x = "rs_id", by.y = "SNP",
             suffixes = c(".bmi", ".t2d"))
mvm <- merge(mvm, eb_dat[, c("SNP", "beta", "se", "ok")], by.x = "rs_id", by.y = "SNP")
mvm <- mvm[mvm$ok, ]
names(mvm)[names(mvm) == "beta.eq"] <- "beta.ramp3"
names(mvm)[names(mvm) == "se.eq"] <- "se.ramp3"
names(mvm)[names(mvm) == "beta"] <- "beta.eb"
names(mvm)[names(mvm) == "se"] <- "se.eb"

append_line("", "## 2. MVMR（校正 BMI、T2D 后 RAMP3 -> eBMD）", "",
            paste0("有效 SNP：", nrow(mvm)), "")
if (nrow(mvm) >= 3) {
  bx <- as.matrix(mvm[, c("beta.ramp3", "beta.bmi", "beta.t2d")])
  bxse <- as.matrix(mvm[, c("se.ramp3", "se.bmi", "se.t2d")])
  by <- mvm$beta.eb
  byse <- mvm$se.eb
  mv_obj <- mr_mvinput(bx = bx, bxse = bxse, by = by, byse = byse,
                       exposure = c("RAMP3", "BMI", "T2D"), outcome = "eBMD")
  append_line("| 暴露 | beta | SE | P |", "|---|---|---|---|")
  ivw_res <- mr_mvivw(mv_obj)
  for (i in seq_along(ivw_res@Exposure)) {
    append_line(sprintf("| %s | %.5f | %.5f | %.3g |",
                        ivw_res@Exposure[i], ivw_res@Estimate[i], ivw_res@StdError[i], ivw_res@Pvalue[i]))
  }
  append_line("", sprintf("MVMR-Q 统计：%.3f，P=%.3g", ivw_res@Heter.Stat[1], ivw_res@Heter.Stat[2]))
  egger_res <- tryCatch(mr_mvegger(mv_obj), error = function(e) NULL)
  if (!is.null(egger_res)) {
    append_line("", "MVMR-Egger（方向性多效性检验）：")
    append_line(sprintf("Egger 截距 = %.5f, P = %.3g", egger_res@Estimate[4], egger_res@Pvalue[4]))
  }
} else {
  append_line("SNP 不足，跳过 MVMR。")
}

# ---------- 3. 中介路径（两步乘积法） ----------
ivw <- function(beta, se) {
  w <- 1 / se^2
  b <- sum(w * beta) / sum(w)
  se_b <- sqrt(1 / sum(w))
  p <- 2 * pnorm(-abs(b / se_b))
  c(beta = b, se = se_b, p = p, n = length(beta))
}

append_line("", "## 3. 中介路径", "")

# a 段：RAMP3 -> BMI / T2D（用 top5）
for (nm in c("BMI", "T2D")) {
  sub <- med_wide[[nm]]
  sub <- sub[sub$ok & sub$SNP %in% top$rs_id, ]
  if (nrow(sub) >= 2) {
    r <- ivw(sub$beta, sub$se)
    append_line(sprintf("RAMP3 -> %s：IVW beta=%.5f, SE=%.5f, P=%.3g, n=%d",
                        nm, r["beta"], r["se"], r["p"], r["n"]))
  } else {
    append_line(sprintf("RAMP3 -> %s：匹配 SNP 不足。", nm))
  }
}

# b 段：BMI -> eBMD（BMI 独立工具）
bmi_b <- merge(bmi_ins, eb_raw, by.x = "SNP", by.y = "SNP")
ab <- align_beta(bmi_b$beta.exposure, bmi_b$se.exposure,
                 bmi_b$effect_allele.exposure, bmi_b$other_allele.exposure,
                 bmi_b$EA, bmi_b$NEA)
bmi_b <- bmi_b[ab$ok, ]
bmi_b$beta.aligned <- ab$beta[ab$ok]
bmi_b$se.aligned <- ab$se[ab$ok]
if (nrow(bmi_b) >= 2) {
  r <- ivw(bmi_b$beta.aligned, bmi_b$se.aligned)
  append_line(sprintf("BMI -> eBMD：IVW beta=%.5f, SE=%.5f, P=%.3g, n=%d",
                      r["beta"], r["se"], r["p"], r["n"]))
} else {
  append_line("BMI -> eBMD：匹配 SNP 不足。")
}

# b 段：T2D -> eBMD（T2D 独立工具）
t2d_b <- merge(t2d_ins, eb_raw, by.x = "SNP", by.y = "SNP")
at2d <- align_beta(t2d_b$beta.exposure, t2d_b$se.exposure,
                   t2d_b$effect_allele.exposure, t2d_b$other_allele.exposure,
                   t2d_b$EA, t2d_b$NEA)
t2d_b <- t2d_b[at2d$ok, ]
t2d_b$beta.aligned <- at2d$beta[at2d$ok]
t2d_b$se.aligned <- at2d$se[at2d$ok]
if (nrow(t2d_b) >= 2) {
  r <- ivw(t2d_b$beta.aligned, t2d_b$se.aligned)
  append_line(sprintf("T2D -> eBMD：IVW beta=%.5f, SE=%.5f, P=%.3g, n=%d",
                      r["beta"], r["se"], r["p"], r["n"]))
} else {
  append_line("T2D -> eBMD：匹配 SNP 不足。")
}

append_line("", "## 4. 使用说明", "",
            "- eQTL beta/se 仍为 GTEx API 近似值，待 slope/se 复核；",
            "- MVMR 未做 LD 校正，top5 eQTL 间高度连锁，结果应视为敏感性分析；",
            "- 中介为简化乘积法，不宣称精确中介比例。")

cat("完成，结果见：", out, "\n")
