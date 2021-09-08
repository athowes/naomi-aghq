#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("replicate_joint-sae-jeff")
# setwd("replicate_joint-sae-jeff")

## devtools::clean_dll("prevartcov")

rstan_options(auto_write = TRUE)
library(parallel)
options(mc.cores = parallel::detectCores())

pars_out <- c(
  "l_rho0",
  "sigma_l_rho",
  "l_alpha0",
  "sigma_l_alpha",
  "rho_i",
  "alpha_i",
  "rho",
  "alpha",
  "bias_rho_i",
  "bias_alpha_i",
  "rmse_rho_i",
  "rmse_alpha_i",
  "bias_l_rho_i",
  "bias_l_alpha_i",
  "rmse_l_rho_i",
  "rmse_l_alpha_i"
)

control_args <- list(adapt_delta = 0.99, max_treedepth = 12)

#' # Test the simulation with 10 simulated datasets

set.seed(5113799)
dat <- replicate(10, sim_data(), FALSE)

sim0 <- mclapply(seq_along(dat), function(i) sampling_summary(object=prevartcov:::stanmodels$model_sim, pars=pars_out, control = control_args, data=c(dat[[i]], model_id=0), refresh=0))
sim1 <- mclapply(seq_along(dat), function(i) sampling_summary(object=prevartcov:::stanmodels$model_sim, pars=pars_out, control = control_args, data=c(dat[[i]], model_id=1), refresh=0))
sim2 <- mclapply(seq_along(dat), function(i) sampling_summary(object=prevartcov:::stanmodels$model_sim, pars=pars_out, control = control_args, data=c(dat[[i]], model_id=2), refresh=0))
sim3 <- mclapply(seq_along(dat), function(i) sampling_summary(object=prevartcov:::stanmodels$model_sim, pars=pars_out, control = control_args, data=c(dat[[i]], model_id=3), refresh=0))
sim4 <- mclapply(seq_along(dat), function(i) sampling_summary(object=prevartcov:::stanmodels$model_sim, pars=pars_out, control = control_args, data=c(dat[[i]], model_id=4), refresh=0))
sim5 <- mclapply(seq_along(dat), function(i) sampling_summary(object=prevartcov:::stanmodels$model_sim, pars=pars_out, control = control_args, data=c(dat[[i]], model_id=5), refresh=0))
sim6 <- mclapply(seq_along(dat), function(i) sampling_summary(object=prevartcov:::stanmodels$model_sim, pars=pars_out, control = control_args, data=c(dat[[i]], model_id=6), refresh=0))

sapply(list(sim0, sim1, sim2, sim3, sim4, sim5, sim6), sapply, "[", "lp__", "Rhat")

rbind(mod0=rowMeans(sapply(sim0, function(x) x[c("bias_rho_i", "bias_alpha_i", "bias_l_rho_i", "bias_l_alpha_i"), 1])),
      mod1=rowMeans(sapply(sim1, function(x) x[c("bias_rho_i", "bias_alpha_i", "bias_l_rho_i", "bias_l_alpha_i"), 1])),
      mod2=rowMeans(sapply(sim2, function(x) x[c("bias_rho_i", "bias_alpha_i", "bias_l_rho_i", "bias_l_alpha_i"), 1])),
      mod3=rowMeans(sapply(sim3, function(x) x[c("bias_rho_i", "bias_alpha_i", "bias_l_rho_i", "bias_l_alpha_i"), 1])),
      mod4=rowMeans(sapply(sim4, function(x) x[c("bias_rho_i", "bias_alpha_i", "bias_l_rho_i", "bias_l_alpha_i"), 1])),
      mod5=rowMeans(sapply(sim5, function(x) x[c("bias_rho_i", "bias_alpha_i", "bias_l_rho_i", "bias_l_alpha_i"), 1])),
      mod6=rowMeans(sapply(sim6, function(x) x[c("bias_rho_i", "bias_alpha_i", "bias_l_rho_i", "bias_l_alpha_i"), 1])))

rbind(mod0=rowMeans(sapply(sim0, function(x) x[c("rmse_rho_i", "rmse_alpha_i", "rmse_l_rho_i", "rmse_l_alpha_i"), 1])),
      mod1=rowMeans(sapply(sim1, function(x) x[c("rmse_rho_i", "rmse_alpha_i", "rmse_l_rho_i", "rmse_l_alpha_i"), 1])),
      mod2=rowMeans(sapply(sim2, function(x) x[c("rmse_rho_i", "rmse_alpha_i", "rmse_l_rho_i", "rmse_l_alpha_i"), 1])),
      mod3=rowMeans(sapply(sim3, function(x) x[c("rmse_rho_i", "rmse_alpha_i", "rmse_l_rho_i", "rmse_l_alpha_i"), 1])),
      mod4=rowMeans(sapply(sim4, function(x) x[c("rmse_rho_i", "rmse_alpha_i", "rmse_l_rho_i", "rmse_l_alpha_i"), 1])),
      mod5=rowMeans(sapply(sim5, function(x) x[c("rmse_rho_i", "rmse_alpha_i", "rmse_l_rho_i", "rmse_l_alpha_i"), 1])),
      mod6=rowMeans(sapply(sim6, function(x) x[c("rmse_rho_i", "rmse_alpha_i", "rmse_l_rho_i", "rmse_l_alpha_i"), 1])))




#' ## Run simulation on DIDE cluster
homedir <- "/Volumes/jwe08"
workdir <- "/Volumes/jeff/prevartcov-simulation"
## dir.create(workdir)

mrc_config <- didehpc::didehpc_config(credentials = "jwe08",
                                      home = homedir,
                                      workdir=workdir,
                                      cluster="mrc",
                                      use_common_lib=FALSE,
                                      use_workers=FALSE)
pkgsrc <- provisionr::package_sources(local="prevartcov_0.0.1.zip")
ctx <- context::context_save(file.path(workdir, "context"), packages="prevartcov", package_sources=pkgsrc)

mrcq <- didehpc::queue_didehpc(ctx, config=mrc_config, initialise=FALSE)


#' Simulate 1000 datasets
#'
set.seed(48470093)
dat <- replicate(1000, sim_data(), FALSE)

#' Do simulation
#'
## mrcq$submit_workers(50)

mrcq$mapply(sampling_summary, data=lapply(dat, c, model_id=0), name="sim0",
            MoreArgs=list(object=prevartcov:::stanmodels$model_sim, iter=5000, control = control_args, pars=pars_out, refresh=0))

mrcq$mapply(sampling_summary, data=lapply(dat, c, model_id=1), name="sim1",
            MoreArgs=list(object=prevartcov:::stanmodels$model_sim, iter=5000, control = control_args, pars=pars_out, refresh=0))
mrcq$mapply(sampling_summary, data=lapply(dat, c, model_id=2), name="sim2",
            MoreArgs=list(object=prevartcov:::stanmodels$model_sim, iter=5000, control = control_args, pars=pars_out, refresh=0))
mrcq$mapply(sampling_summary, data=lapply(dat, c, model_id=3), name="sim3",
            MoreArgs=list(object=prevartcov:::stanmodels$model_sim, iter=5000, control = control_args, pars=pars_out, refresh=0))
mrcq$mapply(sampling_summary, data=lapply(dat, c, model_id=4), name="sim4",
            MoreArgs=list(object=prevartcov:::stanmodels$model_sim, iter=5000, control = control_args, pars=pars_out, refresh=0))
mrcq$mapply(sampling_summary, data=lapply(dat, c, model_id=5), name="sim5",
            MoreArgs=list(object=prevartcov:::stanmodels$model_sim, iter=5000, control = control_args, pars=pars_out, refresh=0))
mrcq$mapply(sampling_summary, data=lapply(dat, c, model_id=6), name="sim6",
            MoreArgs=list(object=prevartcov:::stanmodels$model_sim, iter=5000, control = control_args, pars=pars_out, refresh=0))

## mrcq$stop_workers()

#' Load simulation outputs

sim0 <- mrcq$task_bundle_get("sim0")$results()
sim1 <- mrcq$task_bundle_get("sim1")$results()
sim2 <- mrcq$task_bundle_get("sim2")$results()
sim3 <- mrcq$task_bundle_get("sim3")$results()
sim4 <- mrcq$task_bundle_get("sim4")$results()
sim5 <- mrcq$task_bundle_get("sim5")$results()
sim6 <- mrcq$task_bundle_get("sim6")$results()

simall <- setNames(lapply(paste0("sim", 0:6), get), paste0("mod", 0:6))

#' Check convergence

apply(sapply(simall, sapply, function(x) x["lp__", "Rhat"]), 2, max)
colSums(apply(sapply(simall, sapply, function(x) x["lp__", "Rhat"]), 2, ">", 1.1))
colSums(apply(sapply(simall, sapply, function(x) x["rho", "Rhat"]), 2, ">", 1.1))
colSums(apply(sapply(simall, sapply, function(x) x["alpha", "Rhat"]), 2, ">", 1.1))


#' ## Analysis of simulation outputs

is_covers <- function(sim, data, var, truevar=paste0(var, "_true")){
  est <- sim[grep(paste0("^", var), rownames(sim)),,drop=FALSE]
  est[,"2.5%"] < data[[truevar]] & est[,"97.5%"] > data[[truevar]]
}


#' Estimates for l_prev0 and sigma_l_prev (true value l_prev0 = -2.4 and sigma_l_prev = 0.5)

pars <- c("l_rho0", "sigma_l_rho")

rho0_est <- lapply(simall[-1], sapply, "[", pars, "mean") %>% sapply(rowMeans)

rho0_rmse <- lapply(simall[-1], sapply, "[", pars, "mean") %>%
  sapply(function(x) sqrt(rowMeans((x - c(-2.4, 0.5))^2)))

#' Estimates for l_artcov0 and sigma_l_artcov (true value l_artcov0 = 0.7 and sigma_l_artcov = 0.35)
pars <- c("l_alpha0", "sigma_l_alpha")

alpha0_est <- lapply(simall[4:7], sapply, "[", pars, "mean") %>% sapply(rowMeans)
alpha0_rmse <- lapply(simall[4:7], sapply, "[", pars, "mean") %>%
  sapply(function(x) sqrt(rowMeans((x - c(-2.4, 0.5))^2)))



#' Bias in region-level prevalence and ART coverage
pars <- c("bias_rho_i", "rmse_rho_i", "bias_alpha_i", "rmse_alpha_i")
regest <- lapply(simall[-1], sapply, "[", pars, "mean") %>%
  sapply(rowMeans)

pars <- c("bias_l_rho_i", "rmse_l_rho_i", "bias_l_alpha_i", "rmse_l_alpha_i")
lregest <- round(lapply(simall[-1], sapply, "[", pars, "mean") %>%
                   sapply(rowMeans), 2)

#' # Confidence interval coverage

#' Add true population prevalence to data outputs
dat <- lapply(dat, function(d){d$prev_true <- sum(d$prev_i_true * d$pop_i) / sum(d$pop_i);
d$artcov_true <- sum(d$artcov_i_true * d$prev_i_true * d$pop_i) / sum(d$prev_i_true * d$pop_i); d})
#'

sim_all <- list(sim0=sim0, sim1=sim1, sim2=sim2, sim3=sim3, sim4=sim4, sim5=sim5, sim6=sim6)
prev_i_cov <- sapply(sim_all[-1], function(s) mean(mapply(is_covers, sim=s, data=dat, var="rho_i", truevar = "prev_i_true")))
artcov_i_cov <- sapply(sim_all[-1], function(s) mean(mapply(is_covers, sim=s, data=dat, var="alpha_i", truevar = "artcov_i_true")))


sapply(sim_all[2:7], function(s) mean(mapply(is_covers, sim=s, data=dat, var="l_rho0", truevar = "l_prev0_true")))
sapply(sim_all[4:7], function(s) mean(mapply(is_covers, sim=s, data=dat, var="l_alpha0", truevar="l_artcov0_true")))
sapply(sim_all[2:7], function(s) mean(mapply(is_covers, sim=s, data=dat, var="sigma_l_rho", truevar="sigma_l_prev_true")))
sapply(sim_all[4:7], function(s) mean(mapply(is_covers, sim=s, data=dat, var="sigma_l_alpha", truevar="sigma_l_artcov_true")))


prev_bias <- colMeans(sapply(sim_all[-1], sapply, function(x) x["rho", "mean"]) - sapply(dat, "[[", "prev_true"))
prev_rmse <- sqrt(colMeans((sapply(sim_all[-1], sapply, function(x) x["rho", "mean"]) - sapply(dat, "[[", "prev_true"))^2))

artcov_bias <- colMeans(sapply(sim_all[-1], sapply, function(x) x["alpha", "mean"]) - sapply(dat, "[[", "artcov_true"))
artcov_rmse <- sqrt(colMeans((sapply(sim_all[-1], sapply, function(x) x["alpha", "mean"]) - sapply(dat, "[[", "artcov_true"))^2))

prev_cov <- sapply(sim_all[-1], function(s) mean(mapply(is_covers, sim=s, data=dat, var="rho$", truevar="prev_true")))
artcov_cov <- sapply(sim_all[-1], function(s) mean(mapply(is_covers, sim=s, data=dat, var="alpha$", truevar="artcov_true")))

#' # Make summary table

tab_simout <- rbind(prev_bias, prev_rmse, prev_cov,
                    artcov_bias, artcov_rmse, artcov_cov,
                    regest[1:2,], prev_i_cov,
                    regest[3:4,], artcov_i_cov)

xtable::xtable(tab_simout, digits=1)
