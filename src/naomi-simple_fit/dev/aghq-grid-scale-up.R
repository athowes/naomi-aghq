#' Testing how long it takes to fit the aghq for varying grid sizes

#' (This is the TMB fit)
hypers <- names(fit$par)
length(hypers)

sd_out <- sdreport(obj = fit$obj, par.fixed = fit$par, getJointPrecision = TRUE)

base_grid <- sd_levels_ghe_grid(
  dim = length(hypers),
  level = c(1, 2),
  cut_off = c(0, 1.9),
  sd = sqrt(diag(sd_out$cov.fixed))
)

prod(base_grid$level)

start <- Sys.time()
quad <- fit_aghq(tmb_inputs, basegrid = base_grid)
end <- Sys.time()

end - start

#' TODO change the prod(base_grid$level) upwards
#' Record how long each takes
#' Make graph of relationship
