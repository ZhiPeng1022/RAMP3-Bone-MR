# 拉取 RAMP3 全血全部 eQTL SNP 在血钙/血磷中的关联

suppressMessages(library(TwoSampleMR))

eq <- read.csv("results/summary/data/gtex_ramp3_eqtl.csv")
snps <- unique(eq$snp[eq$tissue == "Whole_Blood"])
cat("SNPs:", length(snps), "\n")

outcomes <- c("ebi-a-GCST90025990", "ebi-a-GCST90025948")
out <- extract_outcome_data(snps = snps, outcomes = outcomes)
out$id.outcome[out$id.outcome == "ebi-a-GCST90025990"] <- "serum_calcium"
out$id.outcome[out$id.outcome == "ebi-a-GCST90025948"] <- "serum_phosphate"

write.csv(out, "results/summary/data/opengwas_outcomes_full.csv", row.names = FALSE)
cat("rows:", nrow(out), "\n")
print(table(out$id.outcome))
