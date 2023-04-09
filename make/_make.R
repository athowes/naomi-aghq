#' prev-anc-art case-study
run_commit_push("prev-anc-art_sim")
run_commit_push("prev-anc-art_model0")
run_commit_push("prev-anc-art_model1")
run_commit_push("prev-anc-art_model2")
run_commit_push("prev-anc-art_model3")
run_commit_push("prev-anc-art_model4")
run_commit_push("prev-anc-art_results")

run_commit_push("epil")
run_commit_push("example_inla-replication")
run_commit_push("example_inla-grid")
run_commit_push("explore_aghq")
run_commit_push("example_naomi")
run_commit_push("explore_posterior-comparison")

#' Naomi model
run_commit_push("naomi")

#' Simplfied Naomi with TMB
id <- orderly::orderly_run("naomi-simple_fit", parameters = list(tmb = TRUE))
orderly::orderly_commit(id)

#' Simplfied Naomi with aghq, and k = 1
id <- orderly::orderly_run("naomi-simple_fit", parameters = list(aghq = TRUE, k = 1, ndConstruction = "product"))
orderly::orderly_commit(id)

#' Simplfied Naomi with aghq, k = 1, and Laplace marginals (a.k.a. adam, for some reason)
id <- orderly::orderly_run("naomi-simple_fit", parameters = list(adam = TRUE))
orderly::orderly_commit(id)

#' Simplfied Naomi with tmbstan, and niter = 1000
#' Note that longer niter are not run locally, and instead use the cluster. See make/hpc_mcmc.R
id <- orderly::orderly_run("naomi-simple_fit", parameters = list(tmbstan = TRUE, hmc_laplace = FALSE))
orderly::orderly_commit(id)

#' Simplfied Naomi with tmbstan, embedded Laplace approximation and niter = 1000
#' Note that longer niter are not run locally, and instead use the cluster. See make/hpc_mcmc.R
id <- orderly::orderly_run("naomi-simple_fit", parameters = list(tmbstan = TRUE, hmc_laplace = TRUE))
orderly::orderly_commit(id)

run_commit_push("naomi-simple_mcmc")
run_commit_push("naomi-simple_compare")
run_commit_push("naomi-simple_model-checks")

#' Statistical methods development

#' Scaling up the grid
run_commit_push("dev_scale-grid")

#' With PCA AGHQ experiments
run_pca_aghq <- function(k, s) {
  id <- orderly::orderly_run("naomi-simple_fit", parameters = list(aghq = TRUE, k = k, s = s, grid_type = "pca"))
  orderly::orderly_commit(id)
}

run_pca_aghq(k = 2, s = 1) #' [x]
run_pca_aghq(k = 2, s = 2) #' [x]
run_pca_aghq(k = 2, s = 3) #' [x]
run_pca_aghq(k = 2, s = 4) #' [x]
run_pca_aghq(k = 2, s = 5) #' [x]
run_pca_aghq(k = 3, s = 1) #' [ ]
run_pca_aghq(k = 3, s = 2) #' [ ]
run_pca_aghq(k = 3, s = 3) #' [ ]
run_pca_aghq(k = 3, s = 4) #' [ ]
run_pca_aghq(k = 3, s = 5) #' [ ]

run_commit_push("naomi-simple_increase-s-k")

#' SINLA development
run_commit_push("dev_sinla") #' Without experiments

#' With experiments
id <- orderly::orderly_run("dev_sinla", param = list(run_experiments = TRUE))
orderly::orderly_commit(id)

#' Documentation and plots
run_commit_push("docs_paper")
run_commit_push("docs_01-04-20-mini")
run_commit_push("docs_01-07-21-stats-epi-group")
run_commit_push("docs_15-11-22-seminar")
run_commit_push("docs_18-04-23-explainer")
run_commit_push("docs_bayescomp-poster")
run_commit_push("plot-tikz_algorithm-flowchart")
run_commit_push("plot-tikz_simplified-naomi")
