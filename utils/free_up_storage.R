#' This orderly stuff is nice, but some of these files are actually getting quite big
#' Who knew?

#' Which of naomi-simple_fit do we need to keep?
orderly::orderly_search("latest(parameter:tmbstan == TRUE && parameter:niter > 50000)", "naomi-simple_fit")
#' "20230414-110613-89cd4244"

#' All of the others should likely be deleted, and start again either without sampling, or with only sampling 1000 rather than 5000
