# 区域级 SuSiE 共定位：eQTLGen RAMP3 核心区域 eQTL vs GEFOS 前臂 BMD
suppressMessages({
  library(susieR)
  library(coloc)
  library(ieugwasr)
})

DATA <- "data/derived"
eq <- read.csv(file.path(DATA, "eqtlgen_core_region_outcomes.csv"), stringsAsFactors = FALSE)
fore <- read.csv(file.path(DATA, "forearm_region_outcomes.csv"), stringsAsFactors = FALSE)
eb <- read.csv(file.path(DATA, "ebmd_chr7_region.csv"), stringsAsFactors = FALSE)

eq <- eq[is.na(eq$proxy.outcome) | !eq$proxy.outcome, ]
fore <- fore[is.na(fore$proxy.outcome) | !fore$proxy.outcome, ]
eq <- eq[!duplicated(eq$SNP), ]
fore <- fore[!duplicated(fore$SNP), ]

m <- merge(eq[, c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome", "other_allele.outcome")],
           fore[, c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome", "other_allele.outcome")],
           by = "SNP", suffixes = c(".eq", ".fore"))
m <- merge(m, eb[, c("RSID", "BP")], by.x = "SNP", by.y = "RSID")

same <- toupper(m$effect_allele.outcome.eq) == toupper(m$effect_allele.outcome.fore) &
  toupper(m$other_allele.outcome.eq) == toupper(m$other_allele.outcome.fore)
flip <- toupper(m$effect_allele.outcome.eq) == toupper(m$other_allele.outcome.fore) &
  toupper(m$other_allele.outcome.eq) == toupper(m$effect_allele.outcome.fore)
m$beta.outcome.fore[flip] <- -m$beta.outcome.fore[flip]
m <- m[same | flip, ]
m <- m[!is.na(m$beta.outcome.eq) & !is.na(m$se.outcome.eq) &
       !is.na(m$beta.outcome.fore) & !is.na(m$se.outcome.fore), ]
m$z_eq <- m$beta.outcome.eq / m$se.outcome.eq
m$z_fore <- m$beta.outcome.fore / m$se.outcome.fore
cat("区域交集 SNP：", nrow(m), "\n")

m$block <- floor(m$BP / 50000)
blocks <- sort(unique(m$block))
cat("LD 块数：", length(blocks), "\n")

ld <- matrix(0, nrow = nrow(m), ncol = nrow(m))
rownames(ld) <- colnames(ld) <- m$SNP
ok_snp <- rep(FALSE, nrow(m))
for (b in blocks) {
  idx <- which(m$block == b)
  if (length(idx) < 2) next
  snps <- m$SNP[idx]
  l <- tryCatch(ieugwasr::ld_matrix(snps, with_alleles = TRUE, pop = "EUR"),
                error = function(e) NULL)
  if (is.null(l)) next
  rs <- sub("_.*", "", rownames(l))
  keep <- rs %in% snps
  l <- l[keep, keep]
  rs <- rs[keep]
  rownames(l) <- colnames(l) <- rs
  ld[rs, rs] <- l
  ok_snp[rownames(ld) %in% rs] <- TRUE
}

m <- m[ok_snp, ]
ld <- ld[rownames(ld) %in% m$SNP, colnames(ld) %in% m$SNP]
cat("有 LD 的 SNP：", nrow(m), "\n")

fit_eq <- susie_rss(m$z_eq, R = ld, L = 10)
fit_fore <- susie_rss(m$z_fore, R = ld, L = 10)
cat("eQTL credible sets：", length(fit_eq$sets$cs), "\n")
cat("forearm credible sets：", length(fit_fore$sets$cs), "\n")

res <- tryCatch(coloc::coloc.susie(fit_eq, fit_fore), error = function(e) NULL)
if (!is.null(res)) {
  print(res$summary)
  write.csv(res$summary, file.path(DATA, "region_susie_coloc_summary.csv"), row.names = FALSE)
}
saveRDS(list(eq_fit = fit_eq, fore_fit = fit_fore, ld = ld, m = m),
        file.path(DATA, "region_susie_coloc_fits.rds"))
cat("完成。\n")
