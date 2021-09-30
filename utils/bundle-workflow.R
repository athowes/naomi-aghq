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
path_bundles <- "bundles"
bundle <- orderly::orderly_bundle_pack(path_bundles, "prev-anc-art_model0-N")

#' A2.
spud <- spud::sharepoint$new("https://imperiallondon.sharepoint.com/")
folder <- spud$folder("HIVInferenceGroup-WP", "Shared Documents/orderly/naomi-inf/bundle-input", verify = TRUE)
folder$upload(path = "bundles/20210930-160113-d37516eb.zip")

#' B: TO-DO
