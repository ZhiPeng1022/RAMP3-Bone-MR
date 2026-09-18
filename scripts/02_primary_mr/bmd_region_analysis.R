# 新 DXA 部位 BMD（GWAS Catalog GCST）区域分析：
# MR（GTEx top5 + eQTLGen）+ 区域 SuSiE 共定位
# 用法：Rscript bmd_region_analysis.R <GCST accession>
suppressMessages({
  library(susieR)
  library(coloc)
  library(ieugwasr)
})

DATA <- "data/derived"
acc <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(acc)) stop("需要 GCST accession")

gc <- read.csv(file.path(DATA, paste0(acc, "_chr7_region.csv")), stringsAsFactors = FALSE)
names(gc)[names(gc) == "beta"] <- "beta.gc"
names(gc)[names(gc) == "se"] <- "se.gc"
names(gc)[names(gc) == "p"] <- "p.gc"
top <- read.delim(file.path(DATA, "ramp3_whole_blood_exposure.txt"), stringsAsFactors = FALSE)
eqgen <- read.csv(file.path(DATA, "eqtlgen_instruments.csv"), stringsAsFactors = FALSE)
eqgen <- eqgen[eqgen$id.exposure == "eqtl-a-ENSG00000122679", ]
eqc <- read.csv(file.path(DATA, "eqtlgen_core_region_outcomes.csv"), stringsAsFactors = FALSE)
eqc <- eqc[is.na(eqc$proxy.outcome) | !eqc$proxy.outcome, ]
eqc <- eqc[!duplicated(eqc$SNP), ]

out_lines <- c(paste0("# ", acc, " 区域分析结果"), "")

ivw_mr <- function(bx, se_x, by, se_y) {
  w <- 1 / se_y^2
  den <- sum(w * bx^2)
  b <- sum(w * bx * by) / den
  se_b <- sqrt(1 / den)
  p <- 2 * pnorm(-abs(b / se_b))
  c(beta = b, se = se_b, p = p, n = length(bx))
}

# ---- MR：GTEx top5 ----
gc_rs <- gc[!is.na(gc$rsid) & gc$rsid != "", ]
m1 <- merge(top, gc_rs, by.x = "rs_id", by.y = "rsid")
if (nrow(m1) >= 2) {
  same <- toupper(m1$EA) == toupper(m1$alt) & toupper(m1$NEA) == toupper(m1$ref)
  flip <- toupper(m1$EA) == toupper(m1$ref) & toupper(m1$NEA) == toupper(m1$alt)
  m1$beta[flip] <- -m1$beta[flip]
  m1 <- m1[same | flip, ]
  m1$beta <- as.numeric(m1$beta); m1$se <- as.numeric(m1$se)
  r <- ivw_mr(m1$beta, m1$se, as.numeric(m1$beta.gc), as.numeric(m1$se.gc))
  out_lines <- c(out_lines, sprintf("GTEx top5 IVW: beta=%.5f, SE=%.5f, P=%.3g, n=%d", r["beta"], r["se"], r["p"], r["n"]))
} else {
  out_lines <- c(out_lines, "GTEx top5: 匹配 SNP 不足")
}

# ---- MR：eQTLGen ----
m2 <- merge(eqgen, gc_rs, by.x = "SNP", by.y = "rsid")
if (nrow(m2) >= 2) {
  same <- toupper(m2$effect_allele.exposure) == toupper(m2$EA) &
    toupper(m2$other_allele.exposure) == toupper(m2$NEA)
  flip <- toupper(m2$effect_allele.exposure) == toupper(m2$NEA) &
    toupper(m2$other_allele.exposure) == toupper(m2$EA)
  m2$beta.gc[flip] <- -m2$beta.gc[flip]
  m2 <- m2[same | flip, ]
  r <- ivw_mr(m2$beta.exposure, m2$se.exposure, as.numeric(m2$beta.gc), as.numeric(m2$se.gc))
  out_lines <- c(out_lines, sprintf("eQTLGen IVW: beta=%.5f, SE=%.5f, P=%.3g, n=%d", r["beta"], r["se"], r["p"], r["n"]))
} else {
  out_lines <- c(out_lines, "eQTLGen: 匹配 SNP 不足")
}

# ---- 区域 SuSiE ----
m3 <- merge(eqc, gc_rs, by.x = "SNP", by.y = "rsid")
same <- toupper(m3$effect_allele.outcome) == toupper(m3$EA) &
  toupper(m3$other_allele.outcome) == toupper(m3$NEA)
flip <- toupper(m3$effect_allele.outcome) == toupper(m3$NEA) &
  toupper(m3$other_allele.outcome) == toupper(m3$EA)
m3$beta.gc[flip] <- -m3$beta.gc[flip]
m3 <- m3[same | flip, ]
m3$z_eq <- m3$beta.outcome / m3$se.outcome
m3$z_gc <- as.numeric(m3$beta.gc) / as.numeric(m3$se.gc)
m3 <- m3[is.finite(m3$z_eq) & is.finite(m3$z_gc), ]
out_lines <- c(out_lines, paste0("区域 SuSiE SNP 数：", nrow(m3)))

if (nrow(m3) >= 100) {
  m3$block <- floor(m3$bp / 50000)
  blocks <- sort(unique(m3$block))
  ld <- matrix(0, nrow = nrow(m3), ncol = nrow(m3))
  rownames(ld) <- colnames(ld) <- m3$SNP
  ok <- rep(FALSE, nrow(m3))
  for (b in blocks) {
    idx <- which(m3$block == b)
    if (length(idx) < 2) next
    l <- tryCatch(ieugwasr::ld_matrix(m3$SNP[idx], with_alleles = TRUE, pop = "EUR"),
                  error = function(e) NULL)
    if (is.null(l)) next
    rs <- sub("_.*", "", rownames(l))
    keep <- rs %in% m3$SNP[idx]
    l <- l[keep, keep]; rs <- rs[keep]
    rownames(l) <- colnames(l) <- rs
    ld[rs, rs] <- l
    ok[rownames(ld) %in% rs] <- TRUE
  }
  m3 <- m3[ok, ]
  ld <- ld[rownames(ld) %in% m3$SNP, colnames(ld) %in% m3$SNP]
  out_lines <- c(out_lines, paste0("有 LD 的 SNP：", nrow(m3)))
  if (nrow(m3) >= 100) {
    fit_eq <- susie_rss(m3$z_eq, R = ld, L = 10)
    fit_gc <- susie_rss(m3$z_gc, R = ld, L = 10)
    out_lines <- c(out_lines, paste0("eQTL credible sets：", length(fit_eq$sets$cs)))
    out_lines <- c(out_lines, paste0("BMD credible sets：", length(fit_gc$sets$cs)))
    res <- tryCatch(coloc::coloc.susie(fit_eq, fit_gc), error = function(e) NULL)
    if (!is.null(res)) {
      write.csv(res$summary, file.path(DATA, paste0(acc, "_region_susie_summary.csv")), row.names = FALSE)
      out_lines <- c(out_lines, "coloc.susie summary 已保存")
    } else {
      out_lines <- c(out_lines, "coloc.susie 无汇总输出")
    }
  }
}

writeLines(out_lines, file.path(DATA, paste0(acc, "_region_analysis.md")))
cat(paste(out_lines, collapse = "\n"), "\n")
