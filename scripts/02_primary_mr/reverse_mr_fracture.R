# C 档-1 反向 MR：DXA 部位 BMD -> 骨折风险
# 暴露：GEFOS 前臂 BMD（ieu-a-977）、股骨颈 BMD（ieu-a-980）
# 结局：FinnGen 前臂骨折（finn-b-ST19_FRACT_FOREA）
suppressMessages(library(TwoSampleMR))

DATA <- "data/derived"

exps <- list(
  "Forearm BMD" = "ieu-a-977",
  "Femoral neck BMD" = "ieu-a-980"
)

result_lines <- character()
for (nm in names(exps)) {
  ins <- tryCatch(
    extract_instruments(exps[[nm]], p1 = 5e-8, clump = TRUE, r2 = 0.01, kb = 10000),
    error = function(e) NULL
  )
  if (is.null(ins) || nrow(ins) == 0) {
    result_lines <- c(result_lines, "", paste0("### ", nm, "：无显著工具（p<5e-8）"), "")
    next
  }
  cat(nm, "工具数：", nrow(ins), "\n")
  out <- extract_outcome_data(snps = ins$SNP, outcomes = "finn-b-ST19_FRACT_FOREA")
  m <- harmonise_data(ins, out)
  m <- m[m$mr_keep, ]
  result_lines <- c(result_lines, "", paste0("### ", nm, " -> Forearm fracture（匹配 ", nrow(m), " SNP）"), "")
  if (nrow(m) >= 2) {
    res <- mr(m, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
    for (i in seq_len(nrow(res))) {
      result_lines <- c(result_lines, sprintf("- %s: beta=%.5f, SE=%.5f, P=%.3g",
                                              res$method[i], res$b[i], res$se[i], res$pval[i]))
    }
    write.csv(res, file.path(DATA, paste0("reverse_mr_", gsub(" ", "_", nm), ".csv")), row.names = FALSE)
  } else {
    result_lines <- c(result_lines, "- 匹配 SNP 不足")
  }
}

out_md <- c(
  "# 反向 MR：DXA 部位 BMD -> 前臂骨折（FinnGen）",
  "",
  "> 暴露方向为 BMD 升高；beta<0 表示 BMD 越高、骨折风险越低",
  result_lines
)
writeLines(out_md, "results/summary/reverse_mr_results.md")
cat(paste(result_lines, collapse = "\n"), "\n")
