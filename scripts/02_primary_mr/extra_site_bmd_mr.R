# 额外部位 BMD 梯度验证：
# Surakka 2020 远端前臂 BMD（ebi-a-GCST90013422, n=21,907）
# Medina-Gomez 2018 total body BMD 45-60 岁（ebi-a-GCST005348, n=18,805）
suppressMessages(library(TwoSampleMR))

DATA <- "data/derived"
top <- read.delim(file.path(DATA, "ramp3_whole_blood_exposure.txt"), stringsAsFactors = FALSE)
eqgen <- read.csv(file.path(DATA, "eqtlgen_instruments.csv"), stringsAsFactors = FALSE)
eqgen <- eqgen[eqgen$id.exposure == "eqtl-a-ENSG00000122679", ]

snps <- unique(c(top$rs_id, eqgen$SNP))
out <- extract_outcome_data(
  snps = snps,
  outcomes = c("ebi-a-GCST90013422", "ebi-a-GCST005348")
)
cat("额外部位匹配行数：", nrow(out), "\n")
write.csv(out, file.path(DATA, "extra_site_bmd_outcomes.csv"), row.names = FALSE)

ivw_mr <- function(beta_x, se_x, beta_y, se_y) {
  w <- 1 / se_y^2
  den <- sum(w * beta_x^2)
  b <- sum(w * beta_x * beta_y) / den
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

result_lines <- character()
for (oid in c("ebi-a-GCST90013422", "ebi-a-GCST005348")) {
  od <- out[out$id.outcome == oid, ]
  trait <- unique(od$outcome)[1]
  result_lines <- c(result_lines, "", paste0("### ", trait, " (", oid, ")"), "")
  m1 <- align(top, od, "effect_allele.outcome", "other_allele.outcome", "alt", "ref")
  if (nrow(m1) >= 2) {
    r <- ivw_mr(m1$beta, m1$se, m1$beta.outcome, m1$se.outcome)
    result_lines <- c(result_lines, sprintf("- GTEx top5: IVW beta=%.5f, SE=%.5f, P=%.3g, n=%d",
                                            r["beta"], r["se"], r["p"], r["n"]))
  } else {
    result_lines <- c(result_lines, "- GTEx top5: 匹配 SNP 不足")
  }
  m2 <- align(eqgen, od, "effect_allele.outcome", "other_allele.outcome",
              "effect_allele.exposure", "other_allele.exposure")
  if (nrow(m2) >= 2) {
    r <- ivw_mr(m2$beta.exposure, m2$se.exposure, m2$beta.outcome, m2$se.outcome)
    result_lines <- c(result_lines, sprintf("- eQTLGen: IVW beta=%.5f, SE=%.5f, P=%.3g, n=%d",
                                            r["beta"], r["se"], r["p"], r["n"]))
  } else {
    result_lines <- c(result_lines, "- eQTLGen: 匹配 SNP 不足")
  }
}

out_md <- c(
  "# 额外部位 BMD 梯度验证结果",
  "",
  "> 暴露 = GTEx v10 top5 RAMP3 eQTL 与 eQTLGen RAMP3 全血 eQTL",
  result_lines
)
writeLines(out_md, "results/summary/extra_site_bmd_results.md")
cat(paste(result_lines, collapse = "\n"), "\n")
