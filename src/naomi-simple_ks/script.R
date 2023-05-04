#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("naomi-simple_ks")
# setwd("src/naomi-simple_ks")

rmarkdown::render("ks.Rmd")

beta_rho <- histogram_and_ecdf("beta_rho", i = 1, return_df = TRUE)
saveRDS(beta_rho$df, "beta_rho.rds")

beta_anc_rho <- histogram_and_ecdf("beta_anc_rho", return_df = TRUE)
saveRDS(beta_anc_rho$df, "beta_anc_rho.rds")
