#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("check_sd-estimation")
# setwd("src/check_sd-estimation")

tmb <- readRDS("depends/tmb.rds")
aghq <- readRDS("depends/aghq.rds")

#' For the parameter ui_anc_alpha_x[18] then TMB appears to estimate the SD as
#' being larger than AGHQ does -- let's look into this

par <- "ui_anc_alpha_x"
i <- 29 #' These are the bad ones: 18, 29, 16, 19

tmb_samples <- tmb$fit$sample[[par]][i, ]
aghq_samples <- aghq$quad$sample[[par]][i, ]

df <- data.frame(
  samples = c(tmb_samples, aghq_samples),
  method = c(rep("TMB", length(tmb_samples)), rep("aghq", length(aghq_samples)))
)

#' Confirming that yes the samples from TMB are more dispersed than those from AGHQ
ggplot(df, aes(x = samples, fill = method)) +
  geom_histogram() +
  facet_grid(~method) +
  scale_fill_manual(values = c("#009E73", "#56B4E9")) +
  guides(fill = "none") +
  labs(title = paste0("Samples for ", par, "[", i, "]"), x = "", y = "") +
  theme_minimal()

#' Now look at the AGHQ modes and Hessians

#' First find the index that this latent field element is stored at
j1 <- which(names(aghq$quad$modesandhessians[["mode"]][[1]]) == par)[i]

aghq_modes <- sapply(aghq$quad$modesandhessians[["mode"]], function(x) x[j1])

tmb_mode <- tmb$fit$mode[[par]][i]

#' The 3^8 AGHQ modes and the 1 TMB mode
data.frame(samples = aghq_modes, method = "AGHQ modes") %>%
  ggplot(aes(x = samples, fill = method)) +
    geom_histogram() +
    geom_vline(xintercept = tmb_mode, color = "#56B4E9") +
    facet_grid(~method) +
    scale_fill_manual(values = c("#009E73")) +
    guides(fill = "none") +
    labs(x = "", y = "") +
    theme_minimal()

TMB::compile("naomi_simple.cpp")
dyn.load(TMB::dynlib("naomi_simple"))

hess <- naomi:::sdreport_joint_precision(tmb$fit$obj, tmb$fit$par.fixed)
cov <- solve(hess)

#' This is the marginal standard deviation
j2 <- which(colnames(cov) == par)[i]
sd_j2 <- sqrt(cov[j2, j2])

(plot <- df %>%
  filter(method == "TMB") %>%
  ggplot(aes(x = samples, fill = method)) +
  geom_histogram(aes(y = after_stat(density)), fill = "#56B4E9", alpha = 0.7) +
  stat_function(fun = dnorm, args = list(mean = tmb_mode, sd = sd_j2)) +
  labs(x = "", y = "") +
  theme_minimal())

#' Expose local_sample_tmb
fit <- tmb$fit
nsample <- 1000
rng_seed <- NULL
random_only <- FALSE
verbose <- TRUE

set.seed(rng_seed)
stopifnot(methods::is(fit, "naomi_fit"))
stopifnot(nsample > 1)
to_tape <- TMB:::isNullPointer(fit$obj$env$ADFun$ptr)
if (to_tape) fit$obj$retape(FALSE)

if (verbose) print("Calculating joint precision")
hess <- naomi:::sdreport_joint_precision(fit$obj, fit$par.fixed)

if (verbose) print("Inverting precision for joint covariance")
cov <- solve(hess)

if (!isSymmetric(cov, tol = sqrt(.Machine$double.eps), check.attributes = FALSE)) {
  stop("cov must be a symmetric matrix")
}

if (verbose) print("Drawing sample")
smp <- mvtnorm::rmvnorm(n = nsample, mean = fit$par.full, sigma = cov, method = "eigen", checkSymmetry = FALSE)

plot +
  geom_histogram(
    data = data.frame(samples = smp[, which(colnames(smp) == par)[i]], method = "Manual TMB"),
    aes(x = samples, y = after_stat(density)), fill = "#E69F00", alpha = 0.7
  )

#' Does the mismatched j correspond to the AGHQ posterior?
sd_j1 <- sqrt(cov[j1, j1])

#' Looks possible!
df %>%
  filter(method == "aghq") %>%
  ggplot(aes(x = samples, fill = method)) +
  geom_histogram(aes(y = after_stat(density)), fill = "#009E73", alpha = 0.7) +
  stat_function(fun = dnorm, args = list(mean = tmb_mode, sd = sd_j1)) +
  labs(x = "", y = "") +
  theme_minimal()

#' Where exactly is the mismatch?
length(colnames(cov)) == length(names(aghq$quad$modesandhessians[["mode"]][[1]]))

#' Oh it's that the AGHQ doesn't have the hyperparameters baked in!
aghq_hess <- aghq$quad$modesandhessians[["H"]][[1]]
aghq_cov <- solve(aghq_hess)
sqrt(aghq_cov[j1, j1])
