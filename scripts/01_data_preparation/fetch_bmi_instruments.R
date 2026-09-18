# 拉取 BMI（ieu-b-40）独立遗传工具

suppressMessages(library(TwoSampleMR))

exp <- extract_instruments(outcomes = "ieu-b-40", p1 = 5e-8)
write.csv(exp, "results/summary/data/bmi_instruments.csv", row.names = FALSE)
cat("BMI instruments:", nrow(exp), "\n")
print(head(exp[, c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure",
                   "other_allele.exposure", "eaf.exposure", "pval.exposure")]))
