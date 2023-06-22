#' A. On local machine:
#' 1. Create bundle
#' 2. Upload to sharepoint somewhere with spud
#'
#' B. On windows VM, or other drive connected to DIDEHPC:
#' 1. Pull bundle
#' 2. Run bundle
#' 3. Upload to sharepoint
#'
#' C. On local machine
#' 1. Pull bundle
#' 2. Import into archive

repo <- "naomi-aghq"
report <- "prev-anc-art_sim"
path_bundles <- "bundles"
param <- NULL

#' A1.
bundle <- orderly::orderly_bundle_pack(path = path_bundles, name = report, parameters = param)

#' A2.
spud <- spud::sharepoint$new("https://imperiallondon.sharepoint.com/")
folder <- spud$folder("HIVInferenceGroup-WP", paste0("Shared Documents/orderly/", repo, "/", path_bundles), verify = TRUE)
folder$upload(path = bundle$path)
recent_bundle <- folder$list() %>% filter(created == max(created))

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
  cluster = "fi--didemrchnb"
  # "fi--dideclusthn"
  # "fi--didemrchnb"
)

# src <- conan::conan_sources(
#   packages = "github::awstringer1/aghq"
# )

src <- NULL

ctx <- context::context_save(
  "context",
  packages = packages,
  package_sources = src
)ß

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

#' C
bundle_output_location <- file.path(root, path_bundles, output_path, t$result()$filename)
orderly::orderly_bundle_import(path = bundle_output_location, root = repo_location)

#' It worked! :party:
