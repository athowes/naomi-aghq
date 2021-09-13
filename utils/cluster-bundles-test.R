#' https://www.vaccineimpact.org/orderly/articles/bundles.html

path <- orderly::orderly_example("prev-anc-art_model0")
path_bundles <- tempfile()

bundle <- orderly::orderly_bundle_pack(path_bundles, "prev-anc-art_model0")

workdir <- tempfile()
res <- orderly::orderly_bundle_run(bundle$path, workdir)

orderly::orderly_bundle_import(res$path)

orderly::orderly_list_archive()
