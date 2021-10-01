#' A. On local machine:
#' 1. Create bundle
#' 2. Upload to sharepoint somewhere with spud
#' B. On windows VM:
#' 1. Pull bundle
#' 2. Run bundle
#' 3. Upload to sharepoint
#' C. On local machine
#' 1. Pull bundle
#' 2. Import into archive

#' A1.
path_bundles <- "bundle-input"
bundle <- orderly::orderly_bundle_pack(path_bundles, "prev-anc-art_model0-N")

#' A2.
spud <- spud::sharepoint$new("https://imperiallondon.sharepoint.com/")
folder <- spud$folder("HIVInferenceGroup-WP", "Shared Documents/orderly/naomi-inf/bundle-input", verify = TRUE)
folder$upload(path = bundle$path)

recent_bundle <- folder$list() %>%
  filter(created == max(created))

#' B1.
folder$download(
  path = recent_bundle$name,
  dest = paste0("bundle-input/", recent_bundle$name)
)

#' B2.
orderly_packages <- yaml::read_yaml(file.path("src/prev-anc-art_model0-N/orderly.yml"))$packages
packages <- list(loaded = c("orderly", orderly_packages))

config <- didehpc::didehpc_config(
  workdir = path_bundles,
  credentials = "ath19",
  cluster = "fi--didemrchnb"
  # "fi--dideclusthn"
  # "fi--didemrchnb"
)

src <- conan::conan_sources("awstringer1/aghq")

ctx <- context::context_save(
  "context",
  packages = packages,
  package_sources = src
)

obj <- didehpc::queue_didehpc(ctx, config = config)

t <- obj$enqueue(orderly::orderly_bundle_run(recent_bundle$name, "bundle-output"))

t$status()
t$result()
