#' Set the working directory to the project root
setwd(rprojroot::find_rstudio_root_file())

#' @param i Index of the artefacts to get
archive_to_docs <- function(report, i = 1) {
  #' Artefacts to be moved
  filenames <- yaml::read_yaml(file.path(paste0("src/", report, "/orderly.yml")))$artefacts[[i]]$data$filenames

  #' Latest version in archive
  latest <- orderly::orderly_latest(report)

  #' Copy files over
  files_from <- paste0("archive/", report, "/", latest, "/", filenames)
  files_to <- paste0("docs/", filenames)

  file.copy(from = files_from, to = files_to, overwrite = TRUE)
}

archive_to_docs("docs_paper")
archive_to_docs("docs_01-04-20-mini")
archive_to_docs("docs_01-07-21-stats-epi-group")
archive_to_docs("docs_15-11-22-seminar")
archive_to_docs("docs_xx-12-22-explainer")
archive_to_docs("epil")
archive_to_docs("prev-anc-art_results")
archive_to_docs("example_inla-replication")
archive_to_docs("example_inla-grid")
archive_to_docs("dev_sinla")
archive_to_docs("dev_sinla", i = 2)
archive_to_docs("explore_aghq")
archive_to_docs("example_naomi")
archive_to_docs("explore_posterior-comparison")
archive_to_docs("docs_bayescomp-poster")
archive_to_docs("plot-tikz_algorithm-flowchart")
archive_to_docs("naomi-simple_mcmc")
archive_to_docs("naomi-simple_compare")
