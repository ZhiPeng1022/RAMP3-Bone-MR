# eQTLGen 补充分析：
# RAMP3 -> eBMD（独立全血 eQTL 复核）
# GLP1R / SLC5A2 -> eBMD（阴性/方法对照）
suppressMessages(library(MendelianRandomization))

DATA <- "data/derived"
ins <- read.csv(file.path(DATA, "eqtlgen_instruments.csv"), stringsAsFactors = FALSE)
eb <- read.csv(file.path(DATA, "ebmd_extract_mvmr.csv"), stringsAsFactors = FALSE)

align_beta <- function(beta, se, ea, oa, out_ea, out_oa) {
  ea <- toupper(ea)
  oa <- toupper(oa)
  out_ea <- toupper(out_ea)
  out_oa <- toupper(out_oa)
  same <- ea == out_ea & oa == out_oa
  flip <- ea == out_oa & oa == out_ea
  same[is.na(same)] <- FALSE
  flip[is.na(flip)] <- FALSE
  beta[flip] <- -beta[flip]
  ok <- same | flip
  list(beta = beta, se = se, ok = ok)
}

out <- "results/summary/eqtlgen_supplement_results.md"
writeLines(c("# eQTLGen 补充分析结果", "",
             "> 数据：eQTLGen Phase I（OpenGWAS eqtl-a-*，全血 n=31,684）", ""), out)
append_line <- function(...) cat(paste0(...), "\n", file = out, append = TRUE)

genes <- list(
  RAMP3 = "eqtl-a-ENSG00000122679",
  GLP1R = "eqtl-a-ENSG00000112164",
  SLC5A2 = "eqtl-a-ENSG00000140675"
)

for (nm in names(genes)) {
  sub <- ins[ins$id.exposure == genes[[nm]], ]
  m <- merge(sub, eb, by.x = "SNP", by.y = "RSID")
  a <- align_beta(m$beta.exposure, m$se.exposure,
                  m$effect_allele.exposure, m$other_allele.exposure,
                  m$EA, m$NEA)
  m <- m[a$ok, ]
  m$beta.exposure <- a$beta[a$ok]
  m$se.exposure <- a$se[a$ok]
  m$beta.eb <- m$BETA
  m$se.eb <- m$SE
  append_line("", paste0("## ", nm, " -> eBMD"), "",
              paste0("eQTLGen 工具数：", nrow(sub), "；eBMD 匹配并对齐：", nrow(m)), "")
  if (nrow(m) >= 3) {
    mr_obj <- mr_input(bx = m$beta.exposure, bxse = m$se.exposure,
                       by = m$beta.eb, byse = m$se.eb)
    ivw <- mr_ivw(mr_obj)
    append_line(sprintf("IVW：beta=%.5f, SE=%.5f, P=%.3g",
                        ivw@Estimate, ivw@StdError, ivw@Pvalue))
    egger <- tryCatch(mr_egger(mr_obj), error = function(e) NULL)
    if (!is.null(egger)) {
      append_line(sprintf("MR-Egger 截距：%.5f, P=%.3g",
                          egger@Intercept, egger@Pvalue.Int))
    }
  } else if (nrow(m) == 1) {
    wr <- m$beta.eb / m$beta.exposure
    wr_se <- m$se.eb / abs(m$beta.exposure)
    p <- 2 * pnorm(-abs(wr / wr_se))
    append_line(sprintf("Wald ratio：beta=%.5f, SE=%.5f, P=%.3g", wr, wr_se, p))
  } else {
    append_line("匹配 SNP 不足，跳过。")
  }
  if (nrow(m) > 0) {
    append_line("", "SNP 明细：")
    for (i in seq_len(nrow(m))) {
      append_line(sprintf("%s: eQTL beta=%.4f, eBMD beta=%.5f",
                          m$SNP[i], m$beta.exposure[i], m$beta.eb[i]))
    }
  }
}

append_line("", "## 解读", "",
            "- RAMP3：eQTLGen 独立全血 eQTL 对 eBMD 的复核（与原 GTEx v10 主分析相互印证）；",
            "- GLP1R：已上市减重靶点对照，观察其对 eBMD 是否同样关联；",
            "- SLC5A2：SGLT2 抑制剂靶点对照，工具数少，结果仅作探索。")

cat("完成：", out, "\n")
