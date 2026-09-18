# Steiger 方向性检验：RAMP3 eQTL -> eBMD / 前臂 BMD
suppressMessages(library(TwoSampleMR))

DATA <- "data/derived"
ins <- read.csv(file.path(DATA, "eqtlgen_instruments.csv"), stringsAsFactors = FALSE)
ins <- ins[ins$id.exposure == "eqtl-a-ENSG00000122679", ]

out <- extract_outcome_data(
  snps = ins$SNP,
  outcomes = c("ebi-a-GCST90014022", "ieu-a-977")
)

res_lines <- character()
for (oid in c("ebi-a-GCST90014022", "ieu-a-977")) {
  od <- out[out$id.outcome == oid, ]
  dat <- harmonise_data(ins, od)
  dat <- dat[dat$mr_keep, ]
  if (nrow(dat) < 2) {
    res_lines <- c(res_lines, paste0(oid, ": SNP 不足"))
    next
  }
  st <- directionality_test(dat)
  print(st)
  res_lines <- c(res_lines, paste0("### ", oid, " | n=", nrow(dat)),
                 paste0("- correct causal direction: ", st$correct_causal_direction),
                 paste0("- steiger p: ", st$steiger_pval),
                 paste0("- snp_r2 exposure/outcome: ", st$snp_r2.exposure, " / ", st$snp_r2.outcome))
}

writeLines(c("# Steiger directionality test", "", res_lines),
           file.path(DATA, "steiger_results.md"))
cat(paste(res_lines, collapse = "\n"), "\n")
