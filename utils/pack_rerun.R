niter <- 1000
nthin <- 10
nchains <- 4
ncores <- 4
report <- "naomi-simple_fit"
param <- list(tmbstan = TRUE, niter = niter, nthin = nthin, nchains = nchains, ncores = ncores)

bundle <- orderly::orderly_bundle_pack(path = "naomi-rerun", name = report, parameters = param)
