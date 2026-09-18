# SNP 级森林图：rs756512（SuSiE 前臂共享信号）与 rs1294873（top5 代表）跨骨结局 Wald ratio
suppressMessages(library(ggplot2))

DATA <- "data/derived"
eq <- read.csv(file.path(DATA, "gtex_ramp3_eqtl.csv"), stringsAsFactors = FALSE)
eq <- eq[eq$tissue == "Whole_Blood" & eq$snp %in% c("rs756512", "rs1294873"), ]
eq$ref <- sub("^.*_([ACGT])_([ACGT])_b38$", "\\1", eq$variant_id)
eq$alt <- sub("^.*_([ACGT])_([ACGT])_b38$", "\\2", eq$variant_id)
eq$z <- qnorm(1 - eq$p / 2)
eq$se <- abs(eq$nes) / eq$z

eb <- read.csv(file.path(DATA, "ebmd_extract_mvmr.csv"), stringsAsFactors = FALSE)
site <- read.csv(file.path(DATA, "site_coloc_outcomes.csv"), stringsAsFactors = FALSE)
gefos <- read.csv(file.path(DATA, "gefos_site_bmd_outcomes.csv"), stringsAsFactors = FALSE)

align_out <- function(beta, se, ea, oa, alt, ref) {
  same <- toupper(ea) == toupper(alt) & toupper(oa) == toupper(ref)
  flip <- toupper(ea) == toupper(ref) & toupper(oa) == toupper(alt)
  beta[flip] <- -beta[flip]
  ok <- same | flip
  list(beta = beta, se = se, ok = ok)
}

rows <- list()
for (i in seq_len(nrow(eq))) {
  snp <- eq$snp[i]
  alt <- eq$alt[i]
  ref <- eq$ref[i]
  # eBMD Morris
  e <- eb[eb$RSID == snp, ]
  if (nrow(e) > 0) {
    a <- align_out(e$BETA, e$SE, e$EA, e$NEA, alt, ref)
    if (any(a$ok)) {
      rows[[length(rows) + 1]] <- data.frame(
        SNP = snp, outcome = "UKB eBMD (Morris 2019)",
        wald = a$beta[a$ok][1] / eq$nes[i],
        wald_se = a$se[a$ok][1] / abs(eq$nes[i]), stringsAsFactors = FALSE)
    }
  }
  # site_coloc: forearm / femoral neck
  for (oid in c("ieu-a-977", "ieu-a-980")) {
    s <- site[site$id.outcome == oid & site$SNP == snp, ]
    if (nrow(s) > 0) {
      a <- align_out(s$beta.outcome, s$se.outcome, s$effect_allele.outcome,
                     s$other_allele.outcome, alt, ref)
      if (any(a$ok)) {
        lab <- ifelse(oid == "ieu-a-977", "Forearm BMD (Zheng 2015)", "Femoral neck BMD (Zheng 2015)")
        rows[[length(rows) + 1]] <- data.frame(
          SNP = snp, outcome = lab,
          wald = a$beta[a$ok][1] / eq$nes[i],
          wald_se = a$se[a$ok][1] / abs(eq$nes[i]),
          stringsAsFactors = FALSE)
      }
    }
  }
  # gefos lumbar
  s2 <- gefos[gefos$id.outcome == "ieu-a-982" & gefos$SNP == snp, ]
  if (nrow(s2) > 0) {
    a <- align_out(s2$beta.outcome, s2$se.outcome, s2$effect_allele.outcome,
                   s2$other_allele.outcome, alt, ref)
    if (any(a$ok)) {
      rows[[length(rows) + 1]] <- data.frame(
        SNP = snp, outcome = "Lumbar spine BMD (Zheng 2015)",
        wald = a$beta[a$ok][1] / eq$nes[i],
        wald_se = a$se[a$ok][1] / abs(eq$nes[i]),
        stringsAsFactors = FALSE)
    }
  }
}

d <- do.call(rbind, rows)
d$ci_low <- d$wald - 1.96 * d$wald_se
d$ci_high <- d$wald + 1.96 * d$wald_se
d$outcome <- factor(d$outcome, levels = rev(unique(d$outcome)))

wrap_label <- function(s, width = 16) {
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

p <- ggplot(d, aes(x = wald, y = outcome, color = SNP)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2, position = position_dodge(width = 0.7)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.7)) +
  scale_color_manual(values = c("rs1294873" = "#8FA3B2", "rs756512" = "#BE8E87")) +
  scale_x_continuous("Wald ratio (per SD of RAMP3 expression)", limits = c(-0.45, 0.10), expand = c(0, 0)) +
  scale_y_discrete(labels = function(x) sapply(x, wrap_label)) +
  labs(title = "Key RAMP3 eQTL SNPs across bone outcomes",
       x = "Wald ratio (per SD of RAMP3 expression)", y = NULL, color = "eQTL SNP") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(hjust = 0, vjust = 0.5, lineheight = 1.0, margin = margin(r = 6)),
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    legend.position = "top",
    legend.justification = "center",
    legend.margin = margin(b = 4),
    plot.margin = margin(t = 10, r = 20, b = 10, l = 20)
  )

ggsave("results/summary/figures/fig6_snp_forest.png", p, width = 7.6, height = 4.6, dpi = 300)
write.csv(d, file.path(DATA, "snp_forest_data.csv"), row.names = FALSE)
print(d)
cat("saved fig7_snp_forest.png\n")
