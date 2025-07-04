# The test version
niter <- 1000
nthin <- 10
nchains <- 4
ncores <- 4
report <- "naomi-simple_fit"
param <- list(tmbstan = TRUE, niter = niter, nthin = nthin, nchains = nchains, ncores = ncores)

bundle <- orderly::orderly_bundle_pack(path = "naomi-rerun", name = report, parameters = param)

# The real version
niter <- 120000
nthin <- 40
nchains <- 18
ncores <- 18
report <- "naomi-simple_fit"
param <- list(tmbstan = TRUE, niter = niter, nthin = nthin, nchains = nchains, ncores = ncores)

bundle <- orderly::orderly_bundle_pack(path = "naomi-rerun", name = report, parameters = param)

recieved_bundle <- "20250619-174833-3548c8b4.zip"

orderly::orderly_bundle_import(path = recieved_bundle, root = NULL)
