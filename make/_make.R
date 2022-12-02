src("make/prev-anc-art.R")
run_commit_push("epil")
run_commit_push("example_inla-replication")
run_commit_push("example_inla-grid")
run_commit_push("explore_aghq")
run_commit_push("example_naomi")
run_commit_push("naomi")

#' Without experiments
run_commit_push("dev_sinla")

#' With experiments
id <- orderly::orderly_run("dev_sinla", param = list(run_experiments = TRUE))
orderly::orderly_commit(id)

run_commit_push("docs_paper")
run_commit_push("docs_01-04-20-mini")
run_commit_push("docs_01-07-21-stats-epi-group")
run_commit_push("docs_15-11-22-seminar")
run_commit_push("docs_xx-12-22-explainer")
run_commit_push("explore_posterior-comparison")
run_commit_push("docs_bayescomp-poster")
run_commit_push("plot-tikz_algorithm-flowchart")
