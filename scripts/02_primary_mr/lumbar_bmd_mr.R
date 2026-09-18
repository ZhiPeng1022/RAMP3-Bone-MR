# 独立骨结局复制：腰椎 BMD（Zheng 2015, ieu-a-982, n=28,498）
# 暴露：GTEx v10 top5 RAMP3 eQTL + eQTLGen RAMP3 工具（两组）
suppressMessages(library(TwoSampleMR))

DATA <- "data/derived"
top <- read.delim(file.path(DATA, "ramp3_whole_blood_exposure.txt"), stringsAsFactors = FALSE)
eqgen <- read.csv(file.path(DATA, "eqtlgen_instruments.csv"), stringsAsFactors = FALSE)
eqgen <- eqgen[eqgen$id.exposure == "eqtl-a-ENSG00000122679", ]

snps <- unique(c(top$rs_id, eqgen$SNP))
out <- extract_outcome_data(snps = snps, outcomes = "ieu-a-982")
cat("腰椎 BMD 匹配行数：", nrow(out), "\n")
write.csv(out, file.path(DATA, "lumbar_bmd_outcomes.csv"), row.names = FALSE)

ivw_mr <- function(beta_x, se_x, beta_y, se_y) {
  w <- 1 / se_y^2
  num <- sum(w * beta_x * beta_y)
  den <- sum(w * beta_x^2)
  b <- num / den
  se_b <- sqrt(1 / den)
  p <- 2 * pnorm(-abs(b / se_b))
  c(beta = b, se = se_b, p = p, n = length(beta_x))
}

align <- function(exposure_df, out_df, ea_col, oa_col, alt_col, ref_col) {
  id_col <- if ("rs_id" %in% names(exposure_df)) "rs_id" else "SNP"
  m <- merge(exposure_df, out_df, by.x = id_col, by.y = "SNP")
  same <- toupper(m[[ea_col]]) == toupper(m[[alt_col]]) & toupper(m[[oa_col]]) == toupper(m[[ref_col]])
  flip <- toupper(m[[ea_col]]) == toupper(m[[ref_col]]) & toupper(m[[oa_col]]) == toupper(m[[alt_col]])
  same[is.na(same)] <- FALSE
  flip[is.na(flip)] <- FALSE
  m$beta.outcome[flip] <- -m$beta.outcome[flip]
  m <- m[same | flip, ]
  m
}

cat("=== GTEx v10 top5 ===")
m1 <- align(top, out, "effect_allele.outcome", "other_allele.outcome", "alt", "ref")
if (nrow(m1) >= 2) {
  r <- ivw_mr(m1$beta, m1$se, m1$beta.outcome, m1$se.outcome)
  cat(sprintf("IVW beta=%.5f SE=%.5f P=%.3g n=%d\n", r["beta"], r["se"], r["p"], r["n"]))
} else {
  cat("匹配 SNP 不足\n")
}

cat("=== eQTLGen RAMP3 ===")
m2 <- align(eqgen, out, "effect_allele.outcome", "other_allele.outcome",
            "effect_allele.exposure", "other_allele.exposure")
if (nrow(m2) >= 2) {
  r <- ivw_mr(m2$beta.exposure, m2$se.exposure, m2$beta.outcome, m2$se.outcome)
  cat(sprintf("IVW beta=%.5f SE=%.5f P=%.3g n=%d\n", r["beta"], r["se"], r["p"], r["n"]))
} else {
  cat("匹配 SNP 不足\n")
}

cat("\nSNP 明细：\n")
print(out[, c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome",
              "other_allele.outcome", "pval.outcome")])
