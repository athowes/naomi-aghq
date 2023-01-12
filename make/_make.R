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
id <- orderly::orderly_run("naomi-simple_fit", parameters = list(aghq = TRUE))
orderly::orderly_commit(id)

#' Simplfied Naomi with tmbstan, and niter = 1000
id <- orderly::orderly_run("naomi-simple_fit", parameters = list(tmbstan = TRUE))
orderly::orderly_commit(id)

#' Statistical methods development
run_commit_push("dev_sinla") #' Without experiments

#' With experiments
id <- orderly::orderly_run("dev_sinla", param = list(run_experiments = TRUE))
orderly::orderly_commit(id)

#' Documentation and plots
run_commit_push("docs_paper")
run_commit_push("docs_01-04-20-mini")
run_commit_push("docs_01-07-21-stats-epi-group")
run_commit_push("docs_15-11-22-seminar")
run_commit_push("docs_xx-12-22-explainer")
run_commit_push("docs_bayescomp-poster")
run_commit_push("plot-tikz_algorithm-flowchart")
