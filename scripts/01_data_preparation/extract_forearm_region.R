# 提取核心区域 SNP 在 GEFOS 前臂 BMD（ieu-a-977）中的效应
suppressMessages(library(TwoSampleMR))

DATA <- "data/derived"
eq <- read.csv(file.path(DATA, "eqtlgen_core_region_outcomes.csv"), stringsAsFactors = FALSE)
direct <- is.na(eq$proxy.outcome) | !eq$proxy.outcome
rsids <- unique(eq$SNP[direct])
cat("核心区域 SNP 数：", length(rsids), "\n")

out <- extract_outcome_data(snps = rsids, outcomes = "ieu-a-977")
cat("前臂 BMD 匹配行数：", nrow(out), "\n")
write.csv(out, file.path(DATA, "forearm_region_outcomes.csv"), row.names = FALSE)
if (nrow(out) > 0) {
  cat("P<0.05 行数：", sum(out$pval.outcome < 0.05, na.rm = TRUE), "\n")
}
