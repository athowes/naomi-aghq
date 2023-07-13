#' 1. Create the scaled PCA grid
#' 2. Visualise the scaled PCA grid marginals
#' 3. Alter fit_aghq to work with the scaled PCA or PCA grids
#' 4. Run fit_aghq with s = 7
#' 5. Check point estimates against PCA-AGHQ in naomi-simple_point-estimates

levels <- c(rep(k, s), rep(1, n_hyper - s))
prod(levels)

pca_scaled_base_grid <- mvQuad::createNIGrid(dim = n_hyper, type = "GHe", level = levels)
