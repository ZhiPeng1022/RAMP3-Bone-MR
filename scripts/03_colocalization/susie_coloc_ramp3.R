# SuSiE 多信号共定位：RAMP3 全血 eQTL（GTEx v10, 105 SNP） vs eBMD
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

eb <- read.csv(file.path(DATA, "ebmd_extract_mvmr.csv"), stringsAsFactors = FALSE)
m <- merge(eq[, c("snp", "variant_id", "alt", "ref", "nes", "se")], eb, by.x = "snp", by.y = "RSID")
same <- toupper(m$EA) == toupper(m$alt) & toupper(m$NEA) == toupper(m$ref)
flip <- toupper(m$EA) == toupper(m$ref) & toupper(m$NEA) == toupper(m$alt)
m$BETA[flip] <- -m$BETA[flip]
m <- m[same | flip, ]
m <- m[!duplicated(m$snp), ]

cat("共享 SNP 数：", nrow(m), "\n")
if (nrow(m) < 10) stop("共享 SNP 太少")

ld <- tryCatch(
  ieugwasr::ld_matrix(m$snp, with_alleles = TRUE, pop = "EUR"),
  error = function(e) {
    cat("LD 获取错误：", conditionMessage(e), "\n")
    NULL
  }
)
if (is.null(ld)) quit(status = 1)
cat("LD 矩阵维度：", dim(ld), "\n")
cat("LD rownames 示例：", paste(head(rownames(ld), 3), collapse = ", "), "\n")

ld_rs <- sub("_.*", "", rownames(ld))
keep <- ld_rs %in% m$snp
ld <- ld[keep, keep]
rownames(ld) <- colnames(ld) <- ld_rs[keep]
common <- intersect(rownames(ld), m$snp)
cat("LD 与 SNP 交集：", length(common), "\n")
m <- m[m$snp %in% common, ]
ld <- ld[common, common]

z_eqtl <- m$nes / m$se
z_ebmd <- m$BETA / m$SE
names(z_eqtl) <- m$snp
names(z_ebmd) <- m$snp

fit_eqtl <- susie_rss(z_eqtl, R = ld, L = 5)
fit_ebmd <- susie_rss(z_ebmd, R = ld, L = 5)

cat("eQTL credible sets：", length(fit_eqtl$sets$cs), "\n")
cat("eBMD credible sets：", length(fit_ebmd$sets$cs), "\n")

res <- tryCatch(
  coloc::coloc.susie(fit_eqtl, fit_ebmd),
  error = function(e) {
    cat("coloc.susie 错误：", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(res)) {
  print(res$summary)
  write.csv(res$summary, file.path(DATA, "susie_coloc_summary.csv"), row.names = FALSE)
}

saveRDS(list(eqtl = fit_eqtl, ebmd = fit_ebmd, ld = ld), file.path(DATA, "susie_coloc_fits.rds"))
cat("完成，结果已保存。\n")
