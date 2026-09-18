# C 档-2 Egger 漏斗图：eQTLGen RAMP3 14 SNP -> 前臂 BMD
suppressMessages(library(ggplot2))

DATA <- "data/derived"
eqgen <- read.csv(file.path(DATA, "eqtlgen_instruments.csv"), stringsAsFactors = FALSE)
eqgen <- eqgen[eqgen$id.exposure == "eqtl-a-ENSG00000122679", ]
site <- read.csv(file.path(DATA, "gefos_site_bmd_outcomes.csv"), stringsAsFactors = FALSE)
fore <- site[site$id.outcome == "ieu-a-977", ]

m <- merge(eqgen, fore, by = "SNP")
same <- toupper(m$effect_allele.exposure) == toupper(m$effect_allele.outcome) &
  toupper(m$other_allele.exposure) == toupper(m$other_allele.outcome)
flip <- toupper(m$effect_allele.exposure) == toupper(m$other_allele.outcome) &
  toupper(m$other_allele.exposure) == toupper(m$effect_allele.outcome)
m$beta.outcome[flip] <- -m$beta.outcome[flip]
m <- m[same | flip, ]

m$ratio <- m$beta.outcome / m$beta.exposure
m$ratio_se <- m$se.outcome / abs(m$beta.exposure)
m$precision <- 1 / m$ratio_se

ivw <- sum((1 / m$ratio_se^2) * m$ratio) / sum(1 / m$ratio_se^2)

p <- ggplot(m, aes(x = ratio, y = precision)) +
  geom_vline(xintercept = ivw, linetype = "dashed", color = "#BE8E87") +
  geom_point(size = 2.5, color = "#8FA3B2") +
  labs(title = "Funnel plot: eQTLGen RAMP3 instruments -> forearm BMD",
       x = "SNP Wald ratio", y = "Precision (1/SE)") +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.margin = margin(t = 10, r = 20, b = 10, l = 20)
  )

ggsave("results/summary/figures/fig7_mr_funnel.png", p, width = 6.4, height = 4.4, dpi = 300)
write.csv(m[, c("SNP", "ratio", "ratio_se")], file.path(DATA, "funnel_data.csv"), row.names = FALSE)
cat("IVW ratio:", ivw, "| SNPs:", nrow(m), "\n")
print(m[, c("SNP", "ratio", "ratio_se")])
