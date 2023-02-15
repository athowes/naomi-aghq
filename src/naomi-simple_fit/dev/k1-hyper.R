#' Debugging issue with hyperparameter for k = 1. Found that this is a feature
#' more than a bug. Need to implement new method for sampling hyperparameters
#' in the k = 1 case that is not based on the grid, otherwise it'll just do EB
#' point estimates

start <- Sys.time()

#' The number of hyperparameters is 24, as compared with 31 for the full model
n_hyper <- 24

if(!(ndConstruction %in% c("product", "sparse"))) {
  warning('ndConstuction must be either "product" or "sparse"')
}

#' Fit AGHQ model
quad <- fit_aghq(tmb_inputs, k = 1)

#' Debug very small standard deviations
summary(quad, method = "correct")

#' Add uncertainty
M <- nsample
M <- 5

# Note that with k = 1, sample_marginal just returns the mode for the hypers
samp <- aghq::sample_marginal(quad, M)
