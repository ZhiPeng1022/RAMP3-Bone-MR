options(repos = c(CRAN = "https://cloud.r-project.org"))

install.packages(c(
  "MendelianRandomization",
  "coloc",
  "susieR",
  "ggplot2",
  "patchwork",
  "remotes"
))

remotes::install_github("MRCIEU/TwoSampleMR")
remotes::install_github("MRCIEU/ieugwasr")
