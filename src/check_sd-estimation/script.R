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
j <- which(names(aghq$quad$modesandhessians[["mode"]][[1]]) == par)[i]

aghq_modes <- sapply(aghq$quad$modesandhessians[["mode"]], function(x) x[j])

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

df %>%
  filter(method == "TMB") %>%
  ggplot(aes(x = samples, fill = method)) +
  geom_histogram(aes(y = after_stat(density)), fill = "#56B4E9") +
  stat_function(fun = dnorm, args = list(mean = tmb_mode, sd = sqrt(cov[j, j]))) +
  labs(x = "", y = "") +
  theme_minimal()
