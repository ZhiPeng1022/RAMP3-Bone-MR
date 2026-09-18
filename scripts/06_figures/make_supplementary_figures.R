# 生成论文图表：森林图、组织表达、单细胞表达、结果表
# 输出：MR/figures/*.png + MR/figures/tables/*.csv

suppressMessages({
  library(ggplot2)
  library(patchwork)
})

DATA <- "data/derived"
FIG <- "figures/final"
TAB <- file.path(FIG, "tables")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

# ---------- 图表主题 ----------
theme_pub <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10, color = "black"),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "none",
    panel.grid.major.y = element_blank()
  )

wrap_label <- function(s, width = 18) {
  if (nchar(s) <= width) return(s)
  words <- strsplit(s, " ")[[1]]
  best_diff <- Inf
  best <- NULL
  for (i in 1:(length(words) - 1)) {
    l1 <- paste(words[1:i], collapse = " ")
    l2 <- paste(words[(i + 1):length(words)], collapse = " ")
    diff <- abs(nchar(l1) - nchar(l2))
    if (diff < best_diff) {
      best_diff <- diff
      best <- c(l1, l2)
    }
  }
  paste(best, collapse = "\n")
}

# ---------- 图 1：组织特异性 Wald ratio 森林图 ----------
tissue <- read.csv(file.path(DATA, "ramp3_tissue_wise_ebmd.csv"))
tissue$lab <- c("Adipose", "Liver", "Muscle", "Thyroid", "Whole blood")[match(tissue$tissue, c(
  "Adipose_Subcutaneous", "Liver", "Muscle_Skeletal", "Thyroid", "Whole_Blood"
))]
tissue$lab <- factor(tissue$lab, levels = rev(tissue$lab))
tissue$ci_low <- tissue$wald_ratio - 1.96 * tissue$wald_se
tissue$ci_high <- tissue$wald_ratio + 1.96 * tissue$wald_se

p1 <- ggplot(tissue, aes(x = wald_ratio, y = lab)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.25, color = "#8FA3B2") +
  geom_point(size = 3, color = "#8FA3B2") +
  scale_x_continuous("Wald ratio (per SD of RAMP3 expression)") +
  scale_y_discrete("Tissue") +
  labs(title = "Tissue-specific RAMP3 eQTL effect on eBMD") +
  theme_pub +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(t = 10, r = 20, b = 10, l = 20)
  )

ggsave(file.path(FIG, "fig1_tissue_mr.png"), p1, width = 6.2, height = 4.2, dpi = 300)

# ---------- 图 2：组织表达 x 遗传信号 ----------
expr <- read.csv(file.path(DATA, "gtex_ramp3_calcr_expression.csv"))
expr <- expr[expr$gene == "RAMP3", ]
expr$lab <- c("Adipose", "Liver", "Muscle", "Thyroid", "Whole blood")[match(expr$tissue, c(
  "Adipose_Subcutaneous", "Liver", "Muscle_Skeletal", "Thyroid", "Whole_Blood"
))]
plot_dat <- merge(expr[, c("lab", "median_expression")], tissue[, c("lab", "p")], by = "lab")
plot_dat$neg_log10p <- -log10(plot_dat$p)

p2a <- ggplot(plot_dat, aes(x = reorder(lab, median_expression), y = median_expression)) +
  geom_col(aes(fill = median_expression), width = 0.65) +
  scale_fill_gradient(low = "#C9D1BB", high = "#7F9372", guide = "none") +
  coord_flip() +
  labs(x = NULL, y = "GTEx median TPM", title = "RAMP3 tissue expression") +
  theme_pub

p2b <- ggplot(plot_dat, aes(x = reorder(lab, neg_log10p), y = neg_log10p)) +
  geom_col(aes(fill = neg_log10p), width = 0.65) +
  scale_fill_gradient(low = "#E5D6D0", high = "#BE8E87", guide = "none") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  coord_flip() +
  labs(x = NULL, y = "-log10(P)", title = "eQTL MR against eBMD") +
  theme_pub

p2 <- p2a + p2b + plot_layout(widths = c(1, 1))
ggsave(file.path(FIG, "fig2_tissue_expression.png"), p2, width = 8.6, height = 3.8, dpi = 300)

# ---------- 图 3：多结局森林图 ----------
outcomes <- data.frame(
  outcome = c("Estimated BMD (SD)", "Any fracture (log OR)",
              "Serum calcium (mmol/L)", "Serum phosphate (mmol/L)"),
  beta = c(0.01232, 0.0023, 0.00233, -0.00224),
  se = c(0.00328, 0.0046, 0.00109, 0.00109),
  p = c(0.000172, 0.607, 0.0328, 0.0405),
  stringsAsFactors = FALSE
)
outcomes$ci_low <- outcomes$beta - 1.96 * outcomes$se
outcomes$ci_high <- outcomes$beta + 1.96 * outcomes$se
outcomes$outcome <- factor(outcomes$outcome, levels = rev(outcomes$outcome))
outcomes$p_lab <- ifelse(outcomes$p < 0.001, "<0.001", sprintf("%.3f", outcomes$p))

p3 <- ggplot(outcomes, aes(x = beta, y = outcome)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.25, color = "#8FA3B2") +
  geom_point(size = 3, color = "#8FA3B2") +
  geom_text(aes(x = ci_high + 0.003, label = paste0("P=", p_lab)), hjust = 0, size = 3.2) +
  scale_x_continuous("Effect size (GWAS units)", limits = c(-0.02, 0.035)) +
  scale_y_discrete("Outcome", labels = function(x) sapply(x, function(s) wrap_label(s, width = 16))) +
  labs(title = "RAMP3 expression proxy and skeletal/mineral outcomes") +
  theme_pub +
  theme(
    axis.text.y = element_text(hjust = 0, vjust = 0.5, lineheight = 1.0, margin = margin(r = 3)),
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(t = 16, r = 36, b = 16, l = 36)
  )

ggsave(file.path(FIG, "fig3_outcomes.png"), p3, width = 7.2, height = 4.4, dpi = 300)

# ---------- 图 4：单细胞与破骨表达 ----------
hpa <- read.csv(file.path(DATA, "hpa_ramp3_singlecell.csv"))
hpa <- head(hpa[order(-hpa$nTPM), ], 12)
hpa$cell_type <- factor(hpa$cell_type, levels = rev(hpa$cell_type))

p4a <- ggplot(hpa, aes(x = nTPM, y = cell_type)) +
  geom_col(fill = "#8E8AA6", width = 0.65) +
  scale_y_discrete(labels = function(x) sapply(x, wrap_label)) +
  labs(x = "Human Protein Atlas nTPM", y = NULL, title = "RAMP3 single-cell expression") +
  theme_pub +
  theme(axis.text.y = element_text(hjust = 0, vjust = 0.5, lineheight = 1.0, margin = margin(r = 6)))

oc <- data.frame(
  gene = c("CTSK", "RAMP3", "RAMP2", "CALCB", "RAMP1", "IAPP", "CALCR"),
  low = c(7.5, 2.7, 1.8, 0.2, 0.1, 0.01, 0.03),
  high = c(43.5, 23.1, 20.4, 2.9, 3.4, 0.14, 0.40),
  stringsAsFactors = FALSE
)
oc$gene <- factor(oc$gene, levels = rev(oc$gene))

p4b <- ggplot(oc, aes(x = gene)) +
  geom_segment(aes(xend = gene, y = low, yend = high), color = "#B8A98C", size = 3) +
  geom_point(aes(y = low), color = "#B8A98C", size = 2.2) +
  geom_point(aes(y = high), color = "#B8A98C", size = 2.2) +
  coord_flip() +
  scale_x_discrete(expand = expansion(mult = c(0.02, 0.06))) +
  labs(x = NULL, y = "Positive cells (%)", title = "Osteoclast-rich tissue (GSE303003)") +
  theme_pub

p4 <- p4a + p4b + plot_layout(widths = c(1.3, 1))
ggsave(file.path(FIG, "fig4_singlecell.png"), p4, width = 10.2, height = 5.6, dpi = 300)

# ---------- 表 1：工具变量表 ----------
eq <- read.delim(file.path(DATA, "ramp3_whole_blood_exposure.txt"))
eb <- read.csv(file.path(DATA, "ebmd_extract_mvmr.csv"))
eb_top <- merge(eq, eb, by.x = "rs_id", by.y = "RSID")
same <- toupper(eb_top$alt) == toupper(eb_top$EA) & toupper(eb_top$ref) == toupper(eb_top$NEA)
flip <- toupper(eb_top$alt) == toupper(eb_top$NEA) & toupper(eb_top$ref) == toupper(eb_top$EA)
eb_top$beta_eb <- ifelse(flip, -eb_top$BETA, eb_top$BETA)
table1 <- data.frame(
  SNP = eb_top$rs_id,
  Effect_allele = eb_top$alt,
  Other_allele = eb_top$ref,
  NES = round(eb_top$beta, 5),
  SE = round(eb_top$se, 5),
  F_stat = round((eb_top$beta / eb_top$se)^2, 1),
  eBMD_beta = round(eb_top$beta_eb, 6),
  eBMD_SE = round(eb_top$SE, 6),
  stringsAsFactors = FALSE
)
write.csv(table1, file.path(TAB, "table1_instruments.csv"), row.names = FALSE)

# ---------- 表 2：多结局 ----------
table2 <- data.frame(
  Outcome = c("eBMD (SD)", "Fracture (log OR)", "Serum calcium (mmol/L)",
              "Serum phosphate (mmol/L)", "BMI (kg/m2)", "T2D (log OR)"),
  Method = c("IVW (top5 eQTL)", "IVW (tissue-wise)", "IVW (whole-blood eQTL)",
             "IVW (whole-blood eQTL)", "IVW (whole-blood eQTL)", "IVW (whole-blood eQTL)"),
  n_SNP = c(5, 4, 6, 6, 20, 60),
  beta = c(0.01232, 0.0023, 0.00233, -0.00224, -0.00166, -0.01623),
  SE = c(0.00328, 0.0046, 0.00109, 0.00109, 0.00047, 0.00121),
  P = c(0.000172, 0.607, 0.0328, 0.0405, 0.00044, 5.5e-41),
  stringsAsFactors = FALSE
)
write.csv(table2, file.path(TAB, "table2_outcomes.csv"), row.names = FALSE)

# ---------- 表 3：MVMR ----------
table3 <- data.frame(
  Exposure = c("RAMP3", "BMI", "T2D"),
  beta = c(-0.00210, 0.99502, -0.15330),
  SE = c(0.00538, 0.35070, 0.07149),
  P = c(0.696, 0.00455, 0.032),
  stringsAsFactors = FALSE
)
attr_note <- data.frame(
  Item = c("n_SNP", "Q", "Q_P", "Egger_intercept", "Egger_P"),
  Value = c("14", "4.180", "0.964", "-0.15683", "0.0286"),
  stringsAsFactors = FALSE
)
write.csv(table3, file.path(TAB, "table3_mvmr.csv"), row.names = FALSE)
write.csv(attr_note, file.path(TAB, "table3_mvmr_notes.csv"), row.names = FALSE)

# ---------- 表 4：中介 ----------
table4 <- data.frame(
  Path = c("RAMP3 -> BMI -> eBMD", "RAMP3 -> T2D -> eBMD"),
  a_beta = c(0.00062, 0.01140),
  a_P = c(0.281, 1.23e-5),
  b_beta = c(-0.00058, -0.00236),
  b_P = c(4.02e-12, 0.00104),
  indirect_beta = c(-3.62e-7, -2.7e-5),
  indirect_SE = c(3.39e-7, 1.03e-5),
  stringsAsFactors = FALSE
)
write.csv(table4, file.path(TAB, "table4_mediation.csv"), row.names = FALSE)

cat("figures saved:", list.files(FIG, pattern = "png$"), "\n")
cat("tables saved:", list.files(TAB), "\n")
