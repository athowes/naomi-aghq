repo <- "naomi-aghq"
report <- "naomi-simple_fit"
path_bundles <- "bundles"

#' k = 3 and s = 9: should take 3 hours
param <- list(aghq = TRUE, k = 3, s = 9, grid_type = "pca", sample = TRUE)

#' k = 3 and s = 10: should take 9 hours
# param <- list(aghq = TRUE, k = 3, s = 9, grid_type = "pca", sample = TRUE)

#' k = 5 and s = 6: should take 1.5 hours
# param <- list(aghq = TRUE, k = 3, s = 9, grid_type = "pca", sample = TRUE)

#' k = 5 and s = 7: should take 7.5 hours
# param <- list(aghq = TRUE, k = 3, s = 9, grid_type = "pca", sample = TRUE)

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
  packages = c("github::athowes/aghq@adam-dev", "github::mrc-ide/naomi")
)

ctx <- context::context_save(
  "context",
  packages = packages,
  package_sources = src
)

obj <- didehpc::queue_didehpc(ctx, config = config)

#' Test that queue works correctly: can also check the version of packages, including aghq, from here
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
t <- queue$task_get("0a5191fb0869a717374a327d65289916")

t$status()
t$log()

#' C
bundle_output_location <- file.path(root, path_bundles, output_path, t$result()$filename)
orderly::orderly_bundle_import(path = bundle_output_location, root = repo_location)
