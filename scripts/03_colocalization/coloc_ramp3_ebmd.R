# 正式 coloc（R coloc 包）：RAMP3 全血 eQTL vs eBMD

library(coloc)

dat <- read.csv("results/summary/data/coloc_input_ramp3_ebmd.csv")

d1 <- list(
  snp = dat$snp,
  beta = dat$beta1,
  varbeta = dat$se1^2,
  type = "quant",
  N = 800,
  MAF = dat$maf
)
d2 <- list(
  snp = dat$snp,
  beta = dat$beta2,
  varbeta = dat$se2^2,
  type = "quant",
  N = 426824,
  MAF = dat$maf
)

res <- coloc.abf(d1, d2)
print(res$summary)
