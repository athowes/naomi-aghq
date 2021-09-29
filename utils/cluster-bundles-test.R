# install.packages("drat")
# drat:::add("mrc-ide")
# install.packages("didehpc")

#' https://www.vaccineimpact.org/orderly/articles/bundles.html

path_bundles <- tempfile()

bundle <- orderly::orderly_bundle_pack(path_bundles, "prev-anc-art_model0-N")

workdir <- tempfile()
res <- orderly::orderly_bundle_run(bundle$path, workdir)

orderly::orderly_bundle_import(res$path)

orderly::orderly_list_archive()

#' https://mrc-ide.github.io/didehpc/articles/orderly.html

path_bundles <- tempfile()

bundle <- orderly::orderly_bundle_pack(path_bundles, "prev-anc-art_model0-N")

#' Obtain the list of packages needed
orderly_packages <- yaml::read_yaml(file.path("src/prev-anc-art_model0-N/orderly.yml"))$packages
packages <- list(loaded = c("orderly", orderly_packages))

my_config <- didehpc::didehpc_config(
  workdir = path_bundles,
  credentials = "ath19",
  cluster = "fi--didemrchnb"
  # "fi--dideclusthn"
  # "fi--didemrchnb"
)

ctx <- context::context_save(root, packages = packages)
