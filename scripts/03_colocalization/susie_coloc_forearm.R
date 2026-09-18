# SuSiE 多信号共定位：RAMP3 eQTL vs 前臂 BMD（ieu-a-977）
suppressMessages({
  library(susieR)
  library(coloc)
  library(ieugwasr)
})

DATA <- "data/derived"
eq <- read.csv(file.path(DATA, "gtex_ramp3_eqtl.csv"), stringsAsFactors = FALSE)
eq <- eq[eq$tissue == "Whole_Blood" & !is.na(eq$nes) & !is.na(eq$p), ]
eq <- eq[!duplicated(eq$snp), ]
pat <- "^.*_([ACGT])_([ACGT])_b38$"
eq$ref <- sub(pat, "\\1", eq$variant_id)
eq$alt <- sub(pat, "\\2", eq$variant_id)
eq$z <- qnorm(1 - eq$p / 2)
eq$se <- abs(eq$nes) / eq$z

site <- read.csv(file.path(DATA, "site_coloc_outcomes.csv"), stringsAsFactors = FALSE)
od <- site[site$id.outcome == "ieu-a-977", ]
od <- od[is.na(od$proxy.outcome) | !od$proxy.outcome, ]
od <- od[!duplicated(od$SNP), ]

m <- merge(eq[, c("snp", "alt", "ref", "nes", "se")], od, by.x = "snp", by.y = "SNP")
same <- toupper(m$effect_allele.outcome) == toupper(m$alt) & toupper(m$other_allele.outcome) == toupper(m$ref)
flip <- toupper(m$effect_allele.outcome) == toupper(m$ref) & toupper(m$other_allele.outcome) == toupper(m$alt)
m$beta.outcome[flip] <- -m$beta.outcome[flip]
m <- m[same | flip, ]
m <- m[!duplicated(m$snp), ]
cat("直接匹配 SNP：", nrow(m), "\n")

ld <- ieugwasr::ld_matrix(m$snp, with_alleles = TRUE, pop = "EUR")
ld_rs <- sub("_.*", "", rownames(ld))
keep <- ld_rs %in% m$snp
ld <- ld[keep, keep]
rownames(ld) <- colnames(ld) <- ld_rs[keep]
m <- m[m$snp %in% rownames(ld), ]
ld <- ld[m$snp, m$snp]
cat("LD 匹配 SNP：", nrow(m), "\n")

z_eqtl <- m$nes / m$se
z_fore <- m$beta.outcome / m$se.outcome
names(z_eqtl) <- names(z_fore) <- m$snp

fit_eqtl <- susie_rss(z_eqtl, R = ld, L = 5)
fit_fore <- susie_rss(z_fore, R = ld, L = 5)
cat("eQTL credible sets：", length(fit_eqtl$sets$cs), "\n")
cat("forearm credible sets：", length(fit_fore$sets$cs), "\n")

res <- tryCatch(coloc::coloc.susie(fit_eqtl, fit_fore), error = function(e) NULL)
if (!is.null(res)) {
  print(res$summary)
  write.csv(res$summary, file.path(DATA, "susie_coloc_forearm_summary.csv"), row.names = FALSE)
} else {
  cat("coloc.susie 未返回结果（可能一侧无信号）\n")
}

saveRDS(list(eqtl = fit_eqtl, forearm = fit_fore, ld = ld), file.path(DATA, "susie_coloc_forearm_fits.rds"))
cat("完成。\n")
