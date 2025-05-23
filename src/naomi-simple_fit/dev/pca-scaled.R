#' [x] 1. Create the scaled PCA grid
#' [x] 2. Visualise the scaled PCA grid marginals
#' [ ] 3. Alter fit_aghq to work with the scaled PCA or PCA grids
#' [ ] 4. Run fit_aghq with s = 7
#' [ ] 5. Check point estimates against PCA-AGHQ in naomi-simple_point-estimates

tmb_input <- tmb_inputs_simple
k <- 3
inner_verbose <- FALSE
progress <- NULL
map <- NULL
DLL <- "naomi_simple"

if (DLL == "naomi_simple") {
  stopifnot(inherits(tmb_input, "naomi_simple_tmb_input"))
} else {
  stopifnot(inherits(tmb_input, "naomi_tmb_input"))
}

obj <- local_make_tmb_obj(tmb_input$data, tmb_input$par_init, calc_outputs = 0L, inner_verbose, progress, map, DLL = DLL)
optresults <- optimize_theta(obj, startingvalue = obj$par, control = default_control_tmb())

#' Can now create whatever grid I'd like using the optresults
pca_grid <- create_pca_grid(m = optresults$mode, C = Matrix::forceSymmetric(solve(optresults$hessian)), s = s, k = k)
nrow(pca_grid$nodes)

S <- length(optresults$mode)
whichfirst <- 1
idxorder <- c(whichfirst, (1:S)[-whichfirst])
m <- optresults$mode[idxorder]
H <- optresults$hessian[idxorder, idxorder]
C <- Matrix::forceSymmetric(solve(H))

scaled_pca_grid <- create_scaled_pca_grid(m, C, s = 7, k = 3)

eigenC <- eigen(C)
lambda_CC <- apply(eigenC$vectors, 2, function(v) as.numeric(v %*% C %*% v / v %*% v))

d <- nrow(C)
var <- diag(C)
hyper_names <- names(m)
hyper_name_starts <- stringr::str_extract(hyper_names, "^[^_]+")
hyper_name_types <- match(hyper_name_starts, unique(hyper_name_starts))
type_means <- tapply(var, hyper_name_types, mean)
mean_var <- type_means[hyper_name_types]
Cs <- diag(1 / sqrt(mean_var)) %*% C %*% diag(1 / sqrt(mean_var))
eigenCs <- eigen(Cs)
lambda_CsC <- apply(eigenCs$vectors, 2, function(v) as.numeric(v %*% C %*% v / v %*% v))

data.frame(x = 1:length(lambda_CsC), Cs = cumsum(lambda_CsC) / sum(lambda_CsC), C = cumsum(lambda_CC) / sum(lambda_CC)) %>%
  pivot_longer(cols = c("Cs", "C"), names_to = "method", values_to = "value") %>%
  ggplot(aes(x = x, y = value, col = as.factor(method))) +
  geom_line() +
  scale_color_brewer(type = "qual") +
  geom_hline(yintercept = 0.9, col = "grey", linetype = "dashed") +
  annotate("text", x = 5, y = 0.875, label = "90% of total variation explained", col = "grey") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "The scaled decomposition explains less variance in the covariance matrix",
    x = "PCA dimensions included", y = "Total variation in C explained", col = "Using the eigendecomposition of"
  ) +
  theme_minimal()

#' Import tmbstan from check_hyper-marginals depends
tmbstan <- readRDS("~/Documents/phd/waterloo/naomi-aghq/src/check_hyper-marginals/depends/tmbstan.rds")

tmbstan_samples <- lapply(
  1:length(hyper_names),
  function(i) data.frame("x" = as.numeric(unlist(tmbstan$mcmc$sample[hyper_names[i]])), par = paste0(hyper_names[i]))
) %>%
  bind_rows()

scaled_pca_nodes <- scaled_pca_grid$nodes %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = "par", values_to = "node")

ggplot(tmbstan_samples, aes(x = x)) +
  geom_histogram(fill = "lightgrey", col = "darkgrey", alpha = 0.5) +
  geom_rug(data = scaled_pca_nodes, aes(x = node), col = "#CC79A7", length = unit(0.1, "npc"), alpha = 0.5) +
  scale_y_continuous(expand = c(0.2, 0.2)) +
  facet_wrap(~ par, scales = "free", ncol = 4) +
  labs(title = "Scaled PCA node positions on hyperparameter marginals from NUTS", x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
