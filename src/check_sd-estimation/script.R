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

pdf("problem-samples.pdf", h = 5, w = 6.25)

ggplot(df, aes(x = samples, fill = method)) +
  geom_histogram() +
  facet_grid(~method) +
  scale_fill_manual(values = c("#009E73", "#E69F00", "#56B4E9")) +
  guides(fill = "none") +
  labs(title = paste0("Samples for ", par, "[", i, "]"), x = "", y = "") +
  theme_minimal()

dev.off()

#' First look at the modes

#' First find the index that this latent field element is stored at
#' I use the notation j1 here because the index is different in different situations
j1 <- which(names(aghq_pca$quad$modesandhessians[["mode"]][[1]]) == par)[i]
aghq_pca_modes <- sapply(aghq_pca$quad$modesandhessians[["mode"]], function(x) x[j1])

tmb_mode <- tmb$fit$mode[[par]][i]
aghq_mode <- aghq$quad$modesandhessians[["mode"]][[1]][j1]

pdf("mode-locations.pdf", h = 5, w = 6.25)

data.frame(samples = aghq_pca_modes, method = "AGHQ modes") %>%
  ggplot(aes(x = samples, fill = method)) +
    geom_histogram() +
    geom_vline(xintercept = tmb_mode, color = "#56B4E9") +
    geom_vline(xintercept = aghq_mode, color = "#009E73") +
    scale_fill_manual(values = c("#E69F00")) +
    guides(fill = "none") +
    labs(title = "The 3^8 AGHQ-PCA modes and the 1 AGHQ and TMB mode", x = "", y = "") +
    theme_minimal()

dev.off()

#' Now look at the Hessians

#' For AGHQ, the Hessians do not have the hyperparameters included
aghq_mode <- aghq$quad$modesandhessians[["mode"]][[1]]
aghq_hess <- aghq$quad$modesandhessians[["H"]][[1]]
dim(aghq_hess) #' 467 x 467

aghq_cov <- solve(aghq_hess)

#' Plot the AGHQ samples and theoretical distribution based on aghq_cov
(plot1 <- ggplot(data.frame(x = aghq$quad$sample[[par]][i, ])) +
    geom_histogram(aes(x = x, y = after_stat(density)), fill = "#009E73", alpha = 0.7) +
    stat_function(fun = dnorm, args = list(mean = aghq_mode[j1], sd = sqrt(aghq_cov[j1, j1]))) +
    labs(x = "", y = "", subtitle = "AGHQ k = 1") +
    theme_minimal())

#' For TMB, the Hessians do have the hyperparameters
tmb_mode <- tmb$fit$par.full
j2 <- which(names(tmb_mode) == par)[i]

dyn.load(TMB::dynlib("naomi_simple"))
tmb_hess <- naomi:::sdreport_joint_precision(tmb$fit$obj, tmb$fit$par.fixed)
dim(tmb_hess) #' 491 x 491

tmb_cov <- solve(tmb_hess)

#' Plot the TMB samples and theoretical distribution based on tmb_cov
(plot2 <- ggplot(data.frame(x = tmb$fit$sample[[par]][i, ])) +
    geom_histogram(aes(x = x, y = after_stat(density)), fill = "#56B4E9", alpha = 0.7) +
    stat_function(fun = dnorm, args = list(mean = tmb_mode[j2], sd = sqrt(tmb_cov[j2, j2]))) +
    labs(x = "", y = "", subtitle = "TMB") +
    theme_minimal())

#' Difference between AGHQ and TMB
plot1 + plot2

#' Look at the two Hessians
image(tmb_hess)
image(aghq_hess)

#' Try excluding the hyperparameters from the TMB Hessian
hyper_indices <- which(colnames(tmb_hess) %in% names(tmb$fit$par))
tmb_hess_reduced <- tmb_hess[-hyper_indices, -hyper_indices]
image(tmb_hess_reduced)

tmb_cov_reduced <- solve(tmb_hess_reduced)

#' Consistent with the hypothesis that these are the same matrix up to numerical error?
diff <- tmb_cov_reduced - aghq_cov

data.frame(
  x = 1:length(diff),
  y = as.vector(diff)
) %>%
  ggplot(aes(x = x, y = y)) +
    geom_point() +
    labs(x = "", y = "tmb_cov_reduced - aghq_cov") +
    theme_minimal()

#' AGHQ sampling algorithm: multinomial a node, then MVN from the Hessian of that node (for the latent field)
#' TMB sampling algorithm: MVN from the Hessian of the joint latent field and hyperparameters
#' If I MVN sample from tmb_cov_reduced then what's the outcome?

#' There is some correlation in here
plot(tmb_cov[hyper_indices[1], ])

nsample <- 1000
tmb_samples_reduced <- mvtnorm::rmvnorm(n = nsample, mean = tmb$fit$par.full[-hyper_indices], sigma = tmb_cov_reduced, method = "eigen", checkSymmetry = FALSE)
j3 <- which(colnames(tmb_samples_reduced) == par)[i]

(plot3 <- ggplot(data.frame(x = tmb_samples_reduced[, j3])) +
    geom_histogram(aes(x = x, y = after_stat(density)), fill = "#E69F00", alpha = 0.7) +
    stat_function(fun = dnorm, args = list(mean = tmb$fit$par.full[-hyper_indices][j3], sd = sqrt(tmb_cov_reduced[j3, j3]))) +
    labs(x = "", y = "", subtitle = "TMB (reduced Hessian)") +
    theme_minimal())

pdf("reduced-hessian-draws.pdf", h = 5, w = 6.25)

plot1 + plot2 + plot3 +
  labs(caption = "Takeaway: AGHQ with k = 1 matches TMB if you use the reduced Hessian")

dev.off()

#' Now, want to look at sdreport$jointPrecision and check this is the same as tmb_hess
#' Answer: yes, up to 10e-9
sdreport <- TMB::sdreport(tmb$fit$obj, par.fixed = tmb$fit$par.fixed, getJointPrecision = TRUE)
max(sdreport$jointPrecision - tmb_hess) < 10e-9

#' Checking what naomi:::sdreport_joint_precision is doing
obj <- tmb$fit$obj
par.fixed <- tmb$fit$par.fixed
hessian.fixed <- NULL
bias.correct <- FALSE
bias.correct.control <- list(sd = FALSE, split = NULL, nsplit = NULL)
ignore.parm.uncertainty <- FALSE
skip.delta.method <- FALSE
