# 提取 eQTLGen RAMP3 核心区域（chr7 b37 45.0-45.4 Mb）eQTL
suppressMessages(library(TwoSampleMR))

DATA <- "data/derived"
eb <- read.csv(file.path(DATA, "ebmd_chr7_region.csv"), stringsAsFactors = FALSE)
eb <- eb[eb$BP >= 45000000 & eb$BP <= 45400000, ]
rsids <- unique(eb$RSID[!is.na(eb$RSID) & eb$RSID != ""])
cat("核心区域 SNP 数：", length(rsids), "\n")

out <- extract_outcome_data(snps = rsids, outcomes = "eqtl-a-ENSG00000122679")
cat("eQTLGen 匹配行数：", nrow(out), "\n")
write.csv(out, file.path(DATA, "eqtlgen_core_region_outcomes.csv"), row.names = FALSE)
if (nrow(out) > 0) {
  direct <- is.na(out$proxy.outcome) | !out$proxy.outcome
  cat("直接行数：", sum(direct), "；P<1e-4：", sum(out$pval.outcome < 1e-4, na.rm = TRUE), "\n")
}
