# 正式 R 版复核 v2：
# 主 MR（全血 top5 eQTL）+ MVMR（全血全部 eQTL 交集）+ 中介路径
# 包：MendelianRandomization
# 输出：formal_analysis_results.md

suppressMessages(library(MendelianRandomization))

DATA <- "data/derived"

eq_top <- read.delim(file.path(DATA, "ramp3_whole_blood_exposure.txt"), stringsAsFactors = FALSE)
eq_all_raw <- read.csv(file.path(DATA, "gtex_ramp3_eqtl.csv"), stringsAsFactors = FALSE)
med <- read.csv(file.path(DATA, "opengwas_mediators.csv"), stringsAsFactors = FALSE)
eb <- read.csv(file.path(DATA, "ebmd_extract_mvmr.csv"), stringsAsFactors = FALSE)
bmi_ins <- read.csv(file.path(DATA, "bmi_instruments.csv"), stringsAsFactors = FALSE)
t2d_ins <- read.csv(file.path(DATA, "t2d_instruments.csv"), stringsAsFactors = FALSE)

cat("top5:", nrow(eq_top), "| eq all:", nrow(eq_all_raw), "| BMI ins:", nrow(bmi_ins),
    "| T2D ins:", nrow(t2d_ins), "| ebmd extracted:", nrow(eb), "\n")

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

# 全血全部 eQTL：从 variant_id 提取 alt/ref，用 nes/p 反推 se
eq_all <- eq_all_raw[eq_all_raw$tissue == "Whole_Blood" & !is.na(eq_all_raw$nes), ]
eq_all <- eq_all[!duplicated(eq_all$snp), ]
pat <- "^.*_([ACGT])_([ACGT])_b38$"
eq_all$ref <- sub(pat, "\\1", eq_all$variant_id)
eq_all$alt <- sub(pat, "\\2", eq_all$variant_id)
eq_all$z <- abs(qnorm(eq_all$p / 2))
eq_all$se <- abs(eq_all$nes) / eq_all$z

# eBMD 原始数据
eb_raw <- data.frame(SNP = eb$RSID, EA = eb$EA, NEA = eb$NEA,
                     beta = eb$BETA, se = eb$SE, stringsAsFactors = FALSE)

out <- "results/summary/formal_analysis_results.md"
writeLines(c("# RAMP3 -> eBMD 正式 R 复核结果", "",
             "> 生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ""), out)
append_line <- function(...) cat(paste0(...), "\n", file = out, append = TRUE)

# ---------- 1. 主 MR：top5 -> eBMD ----------
eb_top_align <- align_beta(eb$BETA, eb$SE, eb$EA, eb$NEA,
                           eq_top$alt[match(eb$RSID, eq_top$rs_id)],
                           eq_top$ref[match(eb$RSID, eq_top$rs_id)])
eb_top_dat <- data.frame(SNP = eb$RSID, beta = eb_top_align$beta, se = eb_top_align$se,
                         ok = eb_top_align$ok, stringsAsFactors = FALSE)
top <- merge(eq_top[, c("rs_id", "beta", "se", "alt", "ref")], eb_top_dat, by.x = "rs_id", by.y = "SNP")
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
}

# ---------- 2. MVMR：全血全部 eQTL ∩ BMI ∩ T2D ----------
get_med_aligned <- function(nm) {
  sub <- med[med$id.outcome == nm & (is.na(med$proxy.outcome) | med$proxy.outcome == FALSE), ]
  sub <- sub[!duplicated(sub$SNP), ]
  if (nrow(sub) > 0) {
    i <- match(sub$SNP, eq_all$snp)
    cat(nm, "example:", sub$SNP[1], sub$effect_allele.outcome[1], sub$other_allele.outcome[1],
        "| eq:", eq_all$alt[i[1]], eq_all$ref[i[1]], "\n")
  }
  a <- align_beta(sub$beta.outcome, sub$se.outcome,
                  sub$effect_allele.outcome, sub$other_allele.outcome,
                  eq_all$alt[match(sub$SNP, eq_all$snp)],
                  eq_all$ref[match(sub$SNP, eq_all$snp)])
  data.frame(SNP = sub$SNP, beta = a$beta, se = a$se, ok = a$ok, stringsAsFactors = FALSE)
}

bmi_w <- get_med_aligned("BMI")
t2d_w <- get_med_aligned("T2D")
cat("bmi direct:", nrow(bmi_w), "ok:", sum(bmi_w$ok), "\n")
cat("t2d direct:", nrow(t2d_w), "ok:", sum(t2d_w$ok), "\n")

eb_all_align <- align_beta(eb$BETA, eb$SE, eb$EA, eb$NEA,
                           eq_all$alt[match(eb$RSID, eq_all$snp)],
                           eq_all$ref[match(eb$RSID, eq_all$snp)])
eb_all_dat <- data.frame(SNP = eb$RSID, beta = eb_all_align$beta, se = eb_all_align$se,
                         ok = eb_all_align$ok, stringsAsFactors = FALSE)
cat("eb all ok:", sum(eb_all_dat$ok), "\n")

eqa <- data.frame(SNP = eq_all$snp, beta.eq = eq_all$nes, se.eq = eq_all$se,
                  F = eq_all$approx_F, stringsAsFactors = FALSE)

bmi_m <- merge(eqa, bmi_w[bmi_w$ok, c("SNP", "beta", "se")], by = "SNP")
names(bmi_m)[names(bmi_m) == "beta"] <- "beta.bmi"
names(bmi_m)[names(bmi_m) == "se"] <- "se.bmi"
t2d_m <- merge(eqa, t2d_w[t2d_w$ok, c("SNP", "beta", "se")], by = "SNP")
names(t2d_m)[names(t2d_m) == "beta"] <- "beta.t2d"
names(t2d_m)[names(t2d_m) == "se"] <- "se.t2d"
eb_m <- merge(eqa, eb_all_dat[eb_all_dat$ok, c("SNP", "beta", "se")], by = "SNP")
names(eb_m)[names(eb_m) == "beta"] <- "beta.eb"
names(eb_m)[names(eb_m) == "se"] <- "se.eb"

mvm <- Reduce(function(a, b) merge(a, b, by = "SNP"), list(bmi_m, t2d_m, eb_m))

run_mvmr <- function(dat, label) {
  append_line("", paste0("### MVMR：", label, "（有效 SNP：", nrow(dat), "）"), "")
  if (nrow(dat) < 3) {
    append_line("SNP 不足，跳过。")
    return(NULL)
  }
  bx <- as.matrix(dat[, c("beta.eq", "beta.bmi", "beta.t2d")])
  bxse <- as.matrix(dat[, c("se.eq", "se.bmi", "se.t2d")])
  by <- dat$beta.eb
  byse <- dat$se.eb
  mv_obj <- mr_mvinput(bx = bx, bxse = bxse, by = by, byse = byse,
                       exposure = c("RAMP3", "BMI", "T2D"), outcome = "eBMD")
  append_line("| 暴露 | beta | SE | P |", "|---|---|---|---|")
  ivw_res <- mr_mvivw(mv_obj)
  for (i in seq_along(ivw_res@Exposure)) {
    append_line(sprintf("| %s | %.5f | %.5f | %.3g |",
                        ivw_res@Exposure[i], ivw_res@Estimate[i], ivw_res@StdError[i], ivw_res@Pvalue[i]))
  }
  append_line("", sprintf("MVMR-Q：%.3f，P=%.3g", ivw_res@Heter.Stat[1], ivw_res@Heter.Stat[2]))
  egger_res <- tryCatch(mr_mvegger(mv_obj), error = function(e) NULL)
  if (!is.null(egger_res)) {
    egger_p <- if ("Pvalue" %in% slotNames(class(egger_res))) egger_res@Pvalue else egger_res@Pvalue.Est
    append_line(sprintf("MVMR-Egger 截距：%.5f，P=%.3g",
                        egger_res@Estimate[length(egger_res@Estimate)],
                        egger_p[length(egger_p)]))
  }
  dat
}

mvm_all <- run_mvmr(mvm, "全部 eQTL 交集")
mvm_f10 <- run_mvmr(mvm[mvm$F >= 10, ], "F>=10 子集")

# 补充 SNP 表
append_line("", "### MVMR 输入 SNP 表（全部交集）", "",
            "| SNP | F(eQTL) | eQTL beta | BMI beta | T2D beta | eBMD beta |",
            "|---|---|---|---|---|---|")
for (i in seq_len(nrow(mvm))) {
  append_line(sprintf("| %s | %.1f | %.4f | %.4f | %.4f | %.5f |",
                      mvm$SNP[i], mvm$F[i], mvm$beta.eq[i], mvm$beta.bmi[i],
                      mvm$beta.t2d[i], mvm$beta.eb[i]))
}

# ---------- 3. 中介路径 ----------
ivw <- function(beta, se) {
  w <- 1 / se^2
  b <- sum(w * beta) / sum(w)
  se_b <- sqrt(1 / sum(w))
  p <- 2 * pnorm(-abs(b / se_b))
  c(beta = b, se = se_b, p = p, n = length(beta))
}

append_line("", "## 3. 中介路径", "")

# a 段：RAMP3 -> BMI / T2D（MVMR 交集 SNP）
for (nm in c("BMI", "T2D")) {
  col_b <- paste0("beta.", tolower(nm))
  col_s <- paste0("se.", tolower(nm))
  if (nrow(mvm) >= 2 && col_b %in% names(mvm)) {
    r <- ivw(mvm[[col_b]], mvm[[col_s]])
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
if (nrow(bmi_b) >= 2) {
  r <- ivw(ab$beta[ab$ok], ab$se[ab$ok])
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
if (nrow(t2d_b) >= 2) {
  r <- ivw(at2d$beta[at2d$ok], at2d$se[at2d$ok])
  append_line(sprintf("T2D -> eBMD：IVW beta=%.5f, SE=%.5f, P=%.3g, n=%d",
                      r["beta"], r["se"], r["p"], r["n"]))
} else {
  append_line("T2D -> eBMD：匹配 SNP 不足。")
}

# 间接效应
append_line("", "### 间接效应（简化乘积法）", "")
a_bmi <- ivw(mvm$beta.bmi, mvm$se.bmi)
b_bmi <- ivw(ab$beta[ab$ok], ab$se[ab$ok])
ind_bmi <- a_bmi["beta"] * b_bmi["beta"]
se_ind_bmi <- sqrt(a_bmi["beta"]^2 * b_bmi["se"]^2 + b_bmi["beta"]^2 * a_bmi["se"]^2)
append_line(sprintf("RAMP3 -> BMI -> eBMD：间接 beta=%.3g，SE=%.3g", ind_bmi, se_ind_bmi))

a_t2d <- ivw(mvm$beta.t2d, mvm$se.t2d)
b_t2d <- ivw(at2d$beta[at2d$ok], at2d$se[at2d$ok])
ind_t2d <- a_t2d["beta"] * b_t2d["beta"]
se_ind_t2d <- sqrt(a_t2d["beta"]^2 * b_t2d["se"]^2 + b_t2d["beta"]^2 * a_t2d["se"]^2)
append_line(sprintf("RAMP3 -> T2D -> eBMD：间接 beta=%.3g，SE=%.3g", ind_t2d, se_ind_t2d))

append_line("", "## 4. 使用说明", "",
            "- eQTL beta/se 为 GTEx API 的 nes/p 近似值，待 slope/se 复核；",
            "- MVMR 仅用全血 eQTL 与 BMI/T2D 直接匹配 SNP，弱工具且未做 LD 校正，视为敏感性分析；",
            "- 中介为简化乘积法，不宣称精确中介比例。")

cat("完成，结果见：", out, "\n")
