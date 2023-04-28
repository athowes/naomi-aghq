#' Fit AGHQ model
quad <- fit_aghq(tmb_inputs_simple, k = k)

if (verbose) print("Sampling from aghq")

M <- 1000
verbose <- TRUE

samp <- aghq::sample_marginal(quad, M)

rownames(samp$samps)

par <- "ui_anc_alpha_x"
i <- 29 #' These are the bad ones: 18, 29, 16, 19

j1 <- which(rownames(samp$samps) == par)[i]

(plot1 <- ggplot(data.frame(x = samp$samps[j1, ])) +
  geom_histogram(aes(x = x, y = after_stat(density)), fill = "#009E73", alpha = 0.7) +
  labs(x = "", y = "", title = "AGHQ k = 1") +
  theme_minimal())

aghq_mode <- quad$modesandhessians[["mode"]][[1]]
aghq_hess <- quad$modesandhessians[["H"]][[1]]
aghq_cov <- solve(aghq_hess)

(plot1 <- plot1 +
  stat_function(fun = dnorm, args = list(mean = aghq_mode[j1], sd = sqrt(aghq_cov[j1, j1]))))

#' Fit TMB model
fit <- local_fit_tmb(tmb_inputs_simple, outer_verbose = TRUE, inner_verbose = FALSE, max_iter = 250, progress = NULL, DLL = "naomi_simple")
fit <- local_sample_tmb(fit, random_only = FALSE, M = M)

(plot2 <- ggplot(data.frame(x = fit$sample[[par]][i, ])) +
    geom_histogram(aes(x = x, y = after_stat(density)), fill = "#56B4E9", alpha = 0.7) +
    labs(x = "", y = "", title = "TMB") +
    theme_minimal())

tmb_mode <- fit$par.full
j2 <- which(names(tmb_mode) == par)[i]

tmb_hess <- naomi:::sdreport_joint_precision(fit$obj, fit$par.fixed)
tmb_cov <- solve(tmb_hess)

(plot2 <- plot2 +
  stat_function(fun = dnorm, args = list(mean = tmb_mode[j2], sd = sqrt(tmb_cov[j2, j2]))))

#' Observe difference between AGHQ and TMB
plot1 + plot2

#' Look at the two Hessians
image(tmb_hess)
image(aghq_hess)

#' Try excluding the hyperparameters from the TMB Hessian
hyper_indices <- which(colnames(tmb_hess) %in% names(fit$par))
tmb_hess_reduced <- tmb_hess[-hyper_indices, -hyper_indices]
image(tmb_hess_reduced)

tmb_cov_reduced <- solve(tmb_hess_reduced)

#' Consistent with the hypothesis that these are the same matrix up to numerical error
sum(tmb_cov_reduced - aghq_cov) / length(tmb_cov_reduced)

#' AGHQ sampling algorithm: multinomial a node, then MVN from the Hessian of that node (for the latent field)
#' TMB sampling algorithm: MVN from the Hessian of the joint latent field and hyperparameters
#' If I MVN sample from tmb_cov_reduced then what's the outcome?

#' There is some correlation in here
plot(tmb_cov[hyper_indices[1], ])

tmb_samples_reduced <- mvtnorm::rmvnorm(n = M, mean = fit$par.full[-hyper_indices], sigma = tmb_cov_reduced, method = "eigen", checkSymmetry = FALSE)
j3 <- which(colnames(tmb_samples_reduced) == par)[i]

(plot3 <- ggplot(data.frame(x = tmb_samples_reduced[, j3])) +
    geom_histogram(aes(x = x, y = after_stat(density)), fill = "#E69F00", alpha = 0.7) +
    labs(x = "", y = "", title = "TMB (reduced Hessian)") +
    theme_minimal())

(plot3 <- plot3 +
    stat_function(fun = dnorm, args = list(mean = fit$par.full[-hyper_indices][j3], sd = sqrt(tmb_cov_reduced[j3, j3]))))

plot1 + plot2 + plot3 +
  labs(caption = "AGHQ with k = 1 matches TMB if you use the reduced Hessian")

