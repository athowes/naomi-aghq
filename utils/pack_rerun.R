niter <- 1000
nthin <- 10
report <- "naomi-simple_fit"
param <- list(tmbstan = TRUE, niter = niter, nthin = nthin)

bundle <- orderly::orderly_bundle_pack(path = "naomi-rerun", name = report, parameters = param)

