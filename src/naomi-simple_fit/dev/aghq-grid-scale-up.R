#' Testing how long it takes to fit the aghq for varying grid sizes

#' Fit TMB model
fit <- local_fit_tmb(tmb_inputs_simple, outer_verbose = TRUE, inner_verbose = FALSE, max_iter = 250, progress = NULL, DLL = "naomi_simple")

hypers <- names(fit$par)
length(hypers)

sd_out <- sdreport(obj = fit$obj, par.fixed = fit$par, getJointPrecision = TRUE)
sd <- sqrt(diag(sd_out$cov.fixed))

#' Some hypers have high SD than others in part due to their scale
#' For example if the hyper is on the logit scale then it will have a higher SD
#' than another hyper on the log scale. Is there a way that I could work around
#' this? Would like to know how much variation in the particular hyper effects
#' variation in outputs

data.frame(par = names(sd), sd = unname(sd)) %>%
  ggplot(aes(x = reorder(par, sd), y = sd)) +
    geom_point() +
    coord_flip() +
    lims(y = c(0, 2)) +
    theme_minimal() +
    labs(x = "", y = "SD from TMB")

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
