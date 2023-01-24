#' Running AGHQ with k = 3 and a sparse grid

repo <- "elgm-inf"
report <- "naomi-simple_fit"
path_bundles <- "bundles"
param <- list(aghq = TRUE, k = 3, ndConstruction = "sparse")
# param <- list(aghq = TRUE) #' For testing

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
  cores = 4
  # "fi--dideclusthn"
  # "fi--didemrchnb"
)

#' aghq, naomi from Github
src <- conan::conan_sources(
  packages = c("github::awstringer1/aghq", "github::mrc-ide/naomi")
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
output_path <- "output"

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
t <- queue$task_get("e7a3af5d7848c44958cdd1e738cc821b")

t$status()
t$log()

#' C
bundle_output_location <- file.path(root, path_bundles, output_path, t$result()$filename)
orderly::orderly_bundle_import(path = bundle_output_location, root = repo_location)
