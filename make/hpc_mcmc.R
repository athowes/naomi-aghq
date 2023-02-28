setting <- 2

#' It was suggested by Seth than 4,000 iterations would be a manageable number to keep
#' This corresponds to four chains of 1000 after thinning. Running locally, it takes
#' 3 hours to generate 4,000 x 4 iterations, but the ESS here is very low, sub 5%. To
#' run 20,000 x 4 iterations on the cluster it should take around 13.5 hours, then
#' we could set nthin = 20 to only keep every 20th iteration, ending up with 1,000 x 4.
if(setting == 1) {
  niter <- 20000
  nthin <- 20
}

#' After trying this, I have found the number of effective samples still to be too low.
#' 20,000 x 4 took around 12 hours. So 100,000 x 4 should take below 3 days. Then we can
#' try thinning by a factor of 20 again to keep only over 40th iteration. If this is too
#' big to practically work with, we can thin it again outside the cluster, but I'd prefer
#' to save too many then thin rather than overthin on the cluster.
if(setting == 2) {
  niter <- 100000
  nthin <- 40
}

repo <- "elgm-inf"
report <- "naomi-simple_fit"
path_bundles <- "bundles"
param <- list(tmbstan = TRUE, niter = niter, nthin = nthin)
# param <- list(tmbstan = TRUE) #' For testing

#' A1.
bundle <- orderly::orderly_bundle_pack(path = path_bundles, name = report, parameters = param)

#' A2.
spud <- spud::sharepoint$new("https://imperiallondon.sharepoint.com/")
folder <- spud$folder("HIVInferenceGroup-WP", paste0("Shared Documents/orderly/", repo, "/", path_bundles), verify = TRUE)
folder$upload(path = bundle$path)
recent_bundle <- dplyr::filter(folder$list(), created == max(created))

#' B
root <- "/Volumes/ath19"
setwd(root)
repo_location <- paste0("~/Documents/waterloo/", repo, "/")

#' B1.
folder$download(
  path = recent_bundle$name,
  dest = file.path(root, path_bundles, recent_bundle$name)
)

#' B2.
orderly_packages <- yaml::read_yaml(
  file.path(paste0(repo_location, "src/", report, "/orderly.yml"))
)$packages

packages <- list(loaded = c("orderly", orderly_packages))

config <- didehpc::didehpc_config(
  workdir = path_bundles,
  credentials = "ath19",
  cluster = "fi--didemrchnb",
  template = "24Core",
  cores = 4,
  # "fi--dideclusthn"
  # "fi--didemrchnb"
)

#' aghq, naomi, TMB from Github
src <- conan::conan_sources(
  packages = c("github::awstringer1/aghq", "github::mrc-ide/naomi", "github::athowes/TMB")
)

ctx <- context::context_save(
  "context",
  packages = packages,
  package_sources = src
)

obj <- didehpc::queue_didehpc(ctx, config = config)

#' Test that queue works correctly
t <- obj$enqueue(sessionInfo())
t$status()
t$result()

#' Run larger job
path <- file.path(recent_bundle$name)
output_path <- file.path("/output")

t <- obj$enqueue(orderly::orderly_bundle_run(
  path = path,
  workdir = output_path
))

t$status()
t$result()

#' Come back to jobs after closing R
ctx_info <- context::context_info("context")
recent_ctx <- dplyr::filter(ctx_info, created == max(created))
ctx <- context::context_load(context::context_read(recent_ctx$id, "context"))
queue <- didehpc::queue_didehpc(ctx)
queue$task_list()
t <- queue$task_get("119c24ff184f8a56c29d27420776b3d8")

t$status()
t$log()

#' C
bundle_output_location <- file.path(root, "output", t$result()$filename)
orderly::orderly_bundle_import(path = bundle_output_location, root = repo_location)
