#' prev-anc-art case-study
run_commit_push("prev-anc-art_sim")
run_commit_push("prev-anc-art_model0")
run_commit_push("prev-anc-art_model1")
run_commit_push("prev-anc-art_model2")
run_commit_push("prev-anc-art_model3")
run_commit_push("prev-anc-art_model4")
run_commit_push("prev-anc-art_results")

run_commit_push("epil")
run_commit_push("astro")
run_commit_push("example_inla-replication")
run_commit_push("example_inla-grid")
run_commit_push("explore_aghq")
run_commit_push("example_naomi")
run_commit_push("explore_posterior-comparison")

#' Naomi model
run_commit_push("naomi")

#' Simplfied Naomi with TMB
id <- orderly::orderly_run("naomi-simple_fit", parameters = list(tmb = TRUE, random_only = TRUE, sample = TRUE))
orderly::orderly_commit(id) #' [x]

#' Simplfied Naomi with aghq, k = 1
id <- orderly::orderly_run("naomi-simple_fit", parameters = list(aghq = TRUE, k = 1, sample = TRUE))
orderly::orderly_commit(id) #' [x]

#' Simplfied Naomi with aghq, k = 1, and Laplace marginals (a.k.a. adam, for some reason)
id <- orderly::orderly_run("naomi-simple_fit", parameters = list(adam = TRUE, sample = TRUE))
orderly::orderly_commit(id) #' [x]

#' Simplfied Naomi with tmbstan, and niter = 1000
#' Note that longer niter are not run locally, and instead use the cluster. See make/hpc_mcmc.R
# id <- orderly::orderly_run("naomi-simple_fit", parameters = list(tmbstan = TRUE, hmc_laplace = FALSE))
# orderly::orderly_commit(id)

#' Simplfied Naomi with tmbstan, embedded Laplace approximation and niter = 1000
# id <- orderly::orderly_run("naomi-simple_fit", parameters = list(tmbstan = TRUE, hmc_laplace = TRUE))
# orderly::orderly_commit(id)

#' Statistical methods development

#' Scaling up the grid
run_commit_push("dev_scale-grid")

#' PCA-AGHQ experiments
run_pca_aghq <- function(k, s, sample = FALSE) {
  id <- orderly::orderly_run("naomi-simple_fit", parameters = list(aghq = TRUE, k = k, s = s, grid_type = "pca", sample = sample))
  orderly::orderly_commit(id)
}

#' Runs for scoping out times and calculating normalising constants
run_pca_aghq(k = 2, s = 1) #' [x]
run_pca_aghq(k = 2, s = 2) #' [x]
run_pca_aghq(k = 2, s = 3) #' [x]
run_pca_aghq(k = 2, s = 4) #' [x]
run_pca_aghq(k = 2, s = 5) #' [x]

run_pca_aghq(k = 3, s = 1) #' [x]
run_pca_aghq(k = 3, s = 2) #' [x]
run_pca_aghq(k = 3, s = 3) #' [x]
run_pca_aghq(k = 3, s = 4) #' [x]
run_pca_aghq(k = 3, s = 5) #' [x]
run_pca_aghq(k = 3, s = 6) #' [x]
run_pca_aghq(k = 3, s = 7) #' [x]
run_pca_aghq(k = 3, s = 8, sample = TRUE) #' [x]

run_pca_aghq(k = 5, s = 1) #' [x]
run_pca_aghq(k = 5, s = 2) #' [x]
run_pca_aghq(k = 5, s = 3) #' [x]
run_pca_aghq(k = 5, s = 4) #' [x]
run_pca_aghq(k = 5, s = 5) #' [x]

run_commit_push("naomi-simple_increase-s-k") #' [x]

#' Scaled PCA-AGHQ experiments
run_scaled_pca_aghq <- function(k, s, sample = FALSE) {
  id <- orderly::orderly_run("naomi-simple_fit", parameters = list(aghq = TRUE, k = k, s = s, grid_type = "scaled_pca", sample = sample))
  orderly::orderly_commit(id)
}

#' Run to first test it works, then a reasonable try at good inference
run_scaled_pca_aghq(k = 3, s = 3) #' [x]
run_scaled_pca_aghq(k = 3, s = 8, sample = TRUE) #' [ ]

#' SINLA development
run_commit_push("dev_sinla") #' Without experiments

#' With experiments
id <- orderly::orderly_run("dev_sinla", param = list(run_experiments = TRUE))
orderly::orderly_commit(id)

#' Results
run_commit_push("naomi-simple_mcmc")
run_commit_push("naomi-simple_ks")
run_commit_push("naomi-simple_exceedance")
run_commit_push("naomi-simple_point-estimates")
run_commit_push("naomi-simple_psis")
run_commit_push("naomi-simple_mmd")
run_commit_push("naomi-simple_model-checks")
run_commit_push("naomi-simple_contraction")

#' Checks
run_commit_push("check_hyper-marginals")
run_commit_push("check_sd-estimation")
run_commit_push("check_tmb-output")
run_commit_push("check_pca-aghq")

#' Documentation and plots
run_commit_push("docs_paper")
run_commit_push("plot-tikz_algorithm-flowchart")
run_commit_push("plot-tikz_simplified-naomi")

#' Presentations
run_commit_push("docs_01-04-20-mini")
run_commit_push("docs_01-07-21-stats-epi-group")
run_commit_push("docs_15-11-22-seminar")
run_commit_push("docs_18-04-23-explainer")
run_commit_push("docs_18-04-23-lab-group")
run_commit_push("docs_19-05-23-turing")
run_commit_push("docs_21-06-23-mlgh")

#' Posters
run_commit_push("docs_bayescomp-poster")
run_commit_push("docs_bioinference-poster")
run_commit_push("docs_statml-poster")
