# The test version
niter <- 1000
nthin <- 10
nchains <- 4
ncores <- 4
report <- "naomi-simple_fit"
param <- list(tmbstan = TRUE, niter = niter, nthin = nthin, nchains = nchains, ncores = ncores)

bundle <- orderly::orderly_bundle_pack(path = "naomi-rerun", name = report, parameters = param)

# The real version
niter <- 80000
nthin <- 10
nchains <- 18
ncores <- 18
report <- "naomi-simple_fit"
param <- list(tmbstan = TRUE, niter = niter, nthin = nthin, nchains = nchains, ncores = ncores)

bundle <- orderly::orderly_bundle_pack(path = "naomi-rerun", name = report, parameters = param)

recieved_bundle <- "20250704-143112-4b78e1d3.zip"

orderly::orderly_bundle_import(path = recieved_bundle, root = NULL)

out <- readRDS("archive/naomi-simple_fit/20250704-143112-4b78e1d3/out.rds")
mcmc <- out$mcmc
s <- rstan::summary(mcmc$stanfit)$summary
s[, "n_eff"] |> summary()
s[, "Rhat"] |> summary()
