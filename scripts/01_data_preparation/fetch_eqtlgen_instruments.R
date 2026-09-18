# 提取 eQTLGen（OpenGWAS eqtl-a-*）的 RAMP3 / GLP1R / SLC5A2 全血 eQTL 工具
# 用途：1) RAMP3 主暴露复核；2) 对照受体工具
suppressMessages(library(TwoSampleMR))

ids <- c(
  "eqtl-a-ENSG00000122679",  # RAMP3
  "eqtl-a-ENSG00000112164",  # GLP1R
  "eqtl-a-ENSG00000140675"   # SLC5A2
)

ins <- extract_instruments(
  outcomes = ids,
  p1 = 1e-4,
  clump = TRUE,
  r2 = 0.01,
  kb = 10000
)

write.csv(ins, "results/summary/data/eqtlgen_instruments.csv", row.names = FALSE)
cat("总工具行数：", nrow(ins), "\n")
if (nrow(ins) > 0) {
  print(table(ins$id.exposure))
  print(head(ins[, c("id.exposure", "SNP", "effect_allele.exposure",
                     "other_allele.exposure", "beta.exposure", "se.exposure",
                     "pval.exposure", "eaf.exposure")]))
  cat("beta 范围：", range(ins$beta.exposure, na.rm = TRUE), "\n")
  cat("se 范围：", range(ins$se.exposure, na.rm = TRUE), "\n")
}
