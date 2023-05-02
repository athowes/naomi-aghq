#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("check_sd-estimation")
# setwd("src/check_sd-estimation")

tmb <- readRDS("depends/tmb.rds")
aghq <- readRDS("depends/aghq.rds")
aghq_pca <- readRDS("depends/aghq-pca.rds")

#' For the parameter ui_anc_alpha_x[29] then TMB appears to estimate the SD as
#' being larger than AGHQ does. This is confusing because AGHQ would seem to take
#' into account more variation thatn TMB. So let's look into this

par <- "ui_anc_alpha_x"
i <- 29 #' These are the bad ones: 18, 29, 16, 19

tmb_samples <- tmb$fit$sample[[par]][i, ]
aghq_samples <- aghq$quad$sample[[par]][i, ]
aghq_pca_samples <- aghq_pca$quad$sample[[par]][i, ]

df <- data.frame(
  samples = c(tmb_samples, aghq_samples, aghq_pca_samples),
  method = c(rep("TMB", length(tmb_samples)), rep("aghq", length(aghq_samples)), rep("aghq-PCA", length(aghq_pca_samples)))
)

#' Confirming that yes the samples from TMB are more dispersed than those from AGHQ
ggplot(df, aes(x = samples, fill = method)) +
  geom_histogram() +
  facet_grid(~method) +
  scale_fill_manual(values = c("#009E73", "#E69F00", "#56B4E9")) +
  guides(fill = "none") +
  labs(title = paste0("Samples for ", par, "[", i, "]"), x = "", y = "") +
  theme_minimal()

#' Now look at modes and Hessians

#' First find the index that this latent field element is stored at
j1 <- which(names(aghq_pca$quad$modesandhessians[["mode"]][[1]]) == par)[i]
aghq_pca_modes <- sapply(aghq_pca$quad$modesandhessians[["mode"]], function(x) x[j1])

tmb_mode <- tmb$fit$mode[[par]][i]
aghq_mode <- aghq$quad$modesandhessians[["mode"]][[1]][j1]

#' The 3^8 AGHQ-PCA modes and the 1 AGHQ and TMB mode
data.frame(samples = aghq_modes, method = "AGHQ modes") %>%
  ggplot(aes(x = samples, fill = method)) +
    geom_histogram() +
    geom_vline(xintercept = tmb_mode, color = "#56B4E9") +
    geom_vline(xintercept = aghq_mode, color = "#009E73") +
    scale_fill_manual(values = c("#E69F00")) +
    guides(fill = "none") +
    labs(x = "", y = "") +
    theme_minimal()

#' Now look at the Hessians

#' AGHQ Hessians do not have the hyperparameters baked in!
samp <- aghq::sample_marginal(aghq$quad, 1000)

(plot1 <- ggplot(data.frame(x = samp$samps[j1, ])) +
    geom_histogram(aes(x = x, y = after_stat(density)), fill = "#009E73", alpha = 0.7) +
    labs(x = "", y = "", title = "AGHQ k = 1") +
    theme_minimal())

aghq_mode <- aghq$quad$modesandhessians[["mode"]][[1]]
aghq_hess <- aghq$quad$modesandhessians[["H"]][[1]]
aghq_cov <- solve(aghq_hess)

#' Plot the AGHQ samples and theoretical distribution
(plot1 <- plot1 +
    stat_function(fun = dnorm, args = list(mean = aghq_mode[j1], sd = sqrt(aghq_cov[j1, j1]))))

#' Want to look at sdreport$jointPrecision
sdreport <- sdreport(tmb$fit$obj, par.fixed = fit$par.fixed, getJointPrecision = TRUE)

(plot2 <- ggplot(data.frame(x = tmb$fit$sample[[par]][i, ])) +
    geom_histogram(aes(x = x, y = after_stat(density)), fill = "#56B4E9", alpha = 0.7) +
    labs(x = "", y = "", title = "TMB") +
    theme_minimal())

tmb_mode <- tmb$fit$par.full
j2 <- which(names(tmb_mode) == par)[i]
tmb_hess <- naomi:::sdreport_joint_precision(tmb$fit$obj, tmb$fit$par.fixed)
tmb_cov <- solve(tmb_hess)

#' Plot the TMB samples and theoretical distribution
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

#' Consistent with the hypothesis that these are the same matrix up to numerical error?
plot(tmb_cov_reduced - aghq_cov)

#' AGHQ sampling algorithm: multinomial a node, then MVN from the Hessian of that node (for the latent field)
#' TMB sampling algorithm: MVN from the Hessian of the joint latent field and hyperparameters
#' If I MVN sample from tmb_cov_reduced then what's the outcome?

#' There is some correlation in here
plot(tmb_cov[hyper_indices[1], ])

tmb_samples_reduced <- mvtnorm::rmvnorm(n = nsample, mean = fit$par.full[-hyper_indices], sigma = tmb_cov_reduced, method = "eigen", checkSymmetry = FALSE)
j3 <- which(colnames(tmb_samples_reduced) == par)[i]

(plot3 <- ggplot(data.frame(x = tmb_samples_reduced[, j3])) +
    geom_histogram(aes(x = x, y = after_stat(density)), fill = "#E69F00", alpha = 0.7) +
    labs(x = "", y = "", title = "TMB (reduced Hessian)") +
    theme_minimal())

(plot3 <- plot3 +
    stat_function(fun = dnorm, args = list(mean = fit$par.full[-hyper_indices][j3], sd = sqrt(tmb_cov_reduced[j3, j3]))))

plot1 + plot2 + plot3 +
  labs(caption = "AGHQ with k = 1 matches TMB if you use the reduced Hessian")

#' Checking what naomi:::sdreport_joint_precision is doing (it's the same as sdreport)
obj <- fit$obj
par.fixed <- fit$par.fixed
hessian.fixed <- NULL
bias.correct <- FALSE
bias.correct.control <- list(sd = FALSE, split = NULL, nsplit = NULL)
ignore.parm.uncertainty <- FALSE
skip.delta.method <- FALSE

# ...
