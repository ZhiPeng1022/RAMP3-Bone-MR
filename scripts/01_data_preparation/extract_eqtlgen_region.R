# 用 eBMD chr7 区域 SNP 列表提取 eQTLGen RAMP3 区域 eQTL 效应
suppressMessages(library(TwoSampleMR))

DATA <- "data/derived"
eb <- read.csv(file.path(DATA, "ebmd_chr7_region.csv"), stringsAsFactors = FALSE)
rsids <- unique(eb$RSID)
rsids <- rsids[!is.na(rsids) & rsids != ""]
cat("区域 SNP 数：", length(rsids), "\n")

out <- extract_outcome_data(
  snps = rsids,
  outcomes = "eqtl-a-ENSG00000122679"
)
cat("eQTLGen 区域匹配行数：", nrow(out), "\n")
write.csv(out, file.path(DATA, "eqtlgen_region_outcomes.csv"), row.names = FALSE)
if (nrow(out) > 0) {
  cat("直接（非代理）行数：", sum(is.na(out$proxy.outcome) | !out$proxy.outcome), "\n")
  cat("P<1e-4 行数：", sum(out$pval.outcome < 1e-4, na.rm = TRUE), "\n")
}
