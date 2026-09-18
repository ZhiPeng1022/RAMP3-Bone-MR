# 拉取 T2D 独立工具变量（Xue 2020，ebi-a-GCST006867）
# 用于：RAMP3 -> T2D -> eBMD 中介路径的 b 段
suppressMessages(library(TwoSampleMR))

t2d_ins <- extract_instruments(
  outcomes = "ebi-a-GCST006867",
  p1 = 5e-8,
  clump = TRUE,
  r2 = 0.01,
  kb = 10000
)
write.csv(t2d_ins, "results/summary/data/t2d_instruments.csv", row.names = FALSE)
cat("T2D instruments:", nrow(t2d_ins), "\n")
if (nrow(t2d_ins) > 0) {
  print(head(t2d_ins[, c("SNP", "effect_allele.exposure", "other_allele.exposure",
                          "beta.exposure", "se.exposure", "pval.exposure")]))
}
