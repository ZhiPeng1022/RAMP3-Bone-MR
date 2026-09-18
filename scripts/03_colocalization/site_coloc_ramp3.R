# DXA 部位 BMD 共定位：
# RAMP3 eQTL（GTEx v10 全血 105 SNP） vs forearm (ieu-a-977) / femoral neck (ieu-a-980)
suppressMessages({
  library(coloc)
  library(TwoSampleMR)
})

DATA <- "data/derived"
eq <- read.csv(file.path(DATA, "gtex_ramp3_eqtl.csv"), stringsAsFactors = FALSE)
eq <- eq[eq$tissue == "Whole_Blood" & !is.na(eq$p), ]
eq <- eq[!duplicated(eq$snp), ]

snps <- unique(eq$snp)
out <- extract_outcome_data(snps = snps, outcomes = c("ieu-a-977", "ieu-a-980"))
write.csv(out, file.path(DATA, "site_coloc_outcomes.csv"), row.names = FALSE)

eqtl_n <- 755  # GTEx v10 Whole_Blood 近似样本量

run_coloc <- function(out_df, outcome_n, label) {
  m <- merge(eq[, c("snp", "p")], out_df, by.x = "snp", by.y = "SNP")
  m <- m[!duplicated(m$snp) & !is.na(m$p) & !is.na(m$pval.outcome) & !is.na(m$eaf.outcome), ]
  cat("=== ", label, " shared SNPs:", nrow(m), "\n")
  if (nrow(m) < 20) return(NULL)
  d1 <- list(snp = m$snp, pvalues = m$p, N = eqtl_n, type = "quant", MAF = m$eaf.outcome)
  d2 <- list(snp = m$snp, pvalues = m$pval.outcome, N = outcome_n, type = "quant", MAF = m$eaf.outcome)
  res <- coloc.abf(d1, d2)
  print(res$summary)
  res
}

res_forearm <- run_coloc(out[out$id.outcome == "ieu-a-977", ], 8143, "Forearm BMD")
res_femoral <- run_coloc(out[out$id.outcome == "ieu-a-980", ], 32735, "Femoral neck BMD")

writeLines(c(
  "# RAMP3 eQTL 与 DXA 部位 BMD 共定位结果",
  "",
  "> coloc.abf；eQTL=GTEx v10 全血（N≈755）；MAF 用结局 GWAS EAF 近似",
  "",
  if (!is.null(res_forearm)) paste0("## Forearm BMD (ieu-a-977, n=8,143)\n\nPP.H4=", res_forearm$summary[6]) else "## Forearm BMD：SNP 不足",
  "",
  if (!is.null(res_femoral)) paste0("## Femoral neck BMD (ieu-a-980, n=32,735)\n\nPP.H4=", res_femoral$summary[6]) else "## Femoral neck BMD：SNP 不足"
), file.path(DATA, "site_coloc_results.md"))
cat("完成，结果保存到 site_coloc_results.md\n")
