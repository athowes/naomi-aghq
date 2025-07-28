# The test version
niter <- 1000
nthin <- 10
nchains <- 4
ncores <- 4
report <- "naomi-simple_fit"
param <- list(tmbstan = TRUE, niter = niter, nthin = nthin, nchains = nchains, ncores = ncores)

bundle <- orderly::orderly_bundle_pack(path = "naomi-rerun", name = report, parameters = param)

# The real version
niter <- 40000
nthin <- 20
nchains <- 18
ncores <- 18
report <- "naomi-simple_fit"
param <- list(tmbstan = TRUE, niter = niter, nthin = nthin, nchains = nchains, ncores = ncores)

bundle <- orderly::orderly_bundle_pack(path = "naomi-rerun", name = report, parameters = param)

recieved_bundle <- "20250723-222404-100e319a"

orderly::orderly_bundle_import(path = paste0(recieved_bundle, ".zip"), root = NULL)

out <- readRDS(file.path("archive", "naomi-simple_fit", recieved_bundle, "out.rds"))
mcmc <- out$mcmc
s <- rstan::summary(mcmc$stanfit)$summary
s[, "n_eff"] |> summary()
s[, "Rhat"] |> summary()
