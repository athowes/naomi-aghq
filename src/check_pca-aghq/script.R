#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("check_pca-aghq")
# setwd("src/check_pca-aghq")

#' Create a function to do the PCA rescaling, which also adapts according to the mean and reweights the nodes:

#' @param m Mean vector
#' @param C Covariance matrix
#' @param s Small grid dimension
#' @param k Number of points per small grid dimension
pca_rescale <- function(m, C, s, k) {
  d <- nrow(C)
  stopifnot(d == length(m))
  eigenC <- eigen(C)
  lambda <- eigenC$values
  Lambda <- diag(lambda)
  E <- eigenC$vectors
  E_s <- E[, 1:s]
  gg_s <- mvQuad::createNIGrid(dim = s, type = "GHe", level = k)
  nodes_out <- t(E_s %*% diag(lambda[1:s]^{0.5}, ncol = s) %*% t(mvQuad::getNodes(gg_s)))
  for(j in 1:d) nodes_out[, j] <- nodes_out[, j] + m[j]
  weights_out <- mvQuad::getWeights(gg_s) * as.numeric(mvQuad::getWeights(mvQuad::createNIGrid(dim = d - s, type = "GHe", level = 1)))

  # Putting things into a mvQuad format manually
  gg <- mvQuad::createNIGrid(dim = d, type = "GHe", level = 1)
  gg$level <- rep(NA, times = d)
  gg$ndConstruction <- "PCA"
  gg$nodes <- nodes_out
  gg$weights <- weights_out
  return(gg)
}

#' Testing with simple examples
#' Here we run some unit tests to make sure, visually and numerically, that the
#' methods are working as expected.

#' Gaussian 2D function: a simple 2D Gaussian example
m <- c(1, 1.5)
C <- matrix(c(2, 1, 1, 1), ncol = 2)

obj <- function(theta) {
  mvtnorm::dmvnorm(theta, mean = m, sigma = C, log = TRUE)
}

ff <- list(
  fn = obj,
  gr = function(theta) numDeriv::grad(obj, theta),
  he = function(theta) numDeriv::hessian(obj, theta)
)

grid <- expand.grid(
  theta1 = seq(-2, 5, length.out = 700),
  theta2 = seq(-2, 5, length.out = 700)
)

ground_truth <- cbind(grid, pdf = exp(obj(grid)))

plot <- ggplot(ground_truth, aes(x = theta1, y = theta2, z = pdf)) +
  geom_contour(col = multi.utils::cbpalette()[1]) +
  coord_fixed(xlim = c(-2, 4.5), ylim = c(-2, 4.5), ratio = 1) +
  labs(x = "", y = "") +
  theme_minimal()

plot

#' Define four AGHQ-PCA grids, with `k` running from 1 to 7:
pca_grid_1 <- pca_rescale(m = m, C = C, s = 1, k = 1)
pca_grid_3 <- pca_rescale(m = m, C = C, s = 1, k = 3)
pca_grid_5 <- pca_rescale(m = m, C = C, s = 1, k = 5)
pca_grid_7 <- pca_rescale(m = m, C = C, s = 1, k = 7)

plot_points <- function(gg) {
  plot +
    geom_point(
      data = mvQuad::getNodes(gg) %>%
        as.data.frame() %>%
        mutate(weights = mvQuad::getWeights(gg)),
      aes(x = V1, y = V2, size = weights),
      alpha = 0.8,
      col = multi.utils::cbpalette()[2],
      inherit.aes = FALSE
    ) +
    scale_size_continuous(range = c(1, 2)) +
    labs(x = "", y = "", size = "Weight") +
    theme_minimal()
}

plot_points(pca_grid_3)

#' We can optimise the objective function using the BFGS algorithm to recover the
#' mode and Hessian at the mode. Of course, for this example, these are known in
#' advance as we are dealing with a Gaussian, but in general this is how they will be computed.

opt_bfgs <- aghq::optimize_theta(ff, m, control = default_control(method = "BFGS"))
opt_bfgs$mode
opt_bfgs$hessian
Matrix::forceSymmetric(solve(opt_bfgs$hessian))

#' Perform regular AGHQ with a product grid (again with a range of values for `k`):
norm_bfgs_1 <- aghq::normalize_logpost(opt_bfgs, k = 1)
norm_bfgs_3 <- aghq::normalize_logpost(opt_bfgs, k = 3)
norm_bfgs_5 <- aghq::normalize_logpost(opt_bfgs, k = 5)
norm_bfgs_7 <- aghq::normalize_logpost(opt_bfgs, k = 7)
norm_bfgs <- list(norm_bfgs_1, norm_bfgs_3, norm_bfgs_5, norm_bfgs_7)

#' For the AGHQ-PCA grids
norm_pca_bfgs_1 <- local_normalize_logpost(opt_bfgs, basegrid = pca_grid_1, adapt = FALSE)
norm_pca_bfgs_3 <- local_normalize_logpost(opt_bfgs, basegrid = pca_grid_3, adapt = FALSE)
norm_pca_bfgs_5 <- local_normalize_logpost(opt_bfgs, basegrid = pca_grid_5, adapt = FALSE)
norm_pca_bfgs_7 <- local_normalize_logpost(opt_bfgs, basegrid = pca_grid_7, adapt = FALSE)
norm_pca_bfgs <- list(norm_pca_bfgs_1, norm_pca_bfgs_3, norm_pca_bfgs_5, norm_pca_bfgs_7)

#' Second way of producing the AGHQ-PCA grids directly using `mvQuad`
pca2_grid_1 <- mvQuad::createNIGrid(dim = 2, type = "GHe", level = c(1, 1))
mvQuad::rescale(pca2_grid_1, m = m, C = C, dec.type = 1)

pca2_grid_3 <- mvQuad::createNIGrid(dim = 2, type = "GHe", level = c(3, 1))
mvQuad::rescale(pca2_grid_3, m = m, C = C, dec.type = 1)

pca2_grid_5 <- mvQuad::createNIGrid(dim = 2, type = "GHe", level = c(5, 1))
mvQuad::rescale(pca2_grid_5, m = m, C = C, dec.type = 1)

pca2_grid_7 <- mvQuad::createNIGrid(dim = 2, type = "GHe", level = c(7, 1))
mvQuad::rescale(pca2_grid_7, m = m, C = C, dec.type = 1)

plot_points(pca2_grid_3)

norm_pca2_bfgs_1 <- local_normalize_logpost(opt_bfgs, basegrid = pca2_grid_1, adapt = FALSE)
norm_pca2_bfgs_3 <- local_normalize_logpost(opt_bfgs, basegrid = pca2_grid_3, adapt = FALSE)
norm_pca2_bfgs_5 <- local_normalize_logpost(opt_bfgs, basegrid = pca2_grid_5, adapt = FALSE)
norm_pca2_bfgs_7 <- local_normalize_logpost(opt_bfgs, basegrid = pca2_grid_7, adapt = FALSE)
norm_pca2_bfgs <- list(norm_pca2_bfgs_1, norm_pca2_bfgs_3, norm_pca2_bfgs_5, norm_pca2_bfgs_7)

#' The density is already normalised, so it has a normalising constant of 1, and
#' a log normalising constant of log(1) = 0
truelognormconst <- 0

results <- data.frame(
  "method" = c("Truth", paste0("aghq, k = ", c(1, 3, 5, 7)), paste0("aghq-pca, k = ", c(1, 3, 5, 7)), paste0("aghq-pca2, k = ", c(1, 3, 5, 7))),
  "lognormconst" = c(truelognormconst, sapply(c(norm_bfgs, norm_pca_bfgs, norm_pca2_bfgs), function(x) x$lognormconst)),
  "type" = c("Truth", rep("AGHQ", 4), rep("AGHQ-PCA", 4), rep("AGHQ-PCA2", 4))
)

ggplot(results, aes(x = method, y = lognormconst, col = as.factor(type))) +
  geom_point(size = 1) +
  labs(x = "Method", y = "Log normalising constant", col = "") +
  scale_color_manual(values = multi.utils::cbpalette()) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 70, vjust = 0.5, hjust = 0.5))

#' Non-Gaussian 2D function
#' A more complicated two dimensional function from Alex's unit tests
#' https://github.com/awstringer1/aghq/blob/master/tests/testthat/setup-01-optimization.R
set.seed(84343124)

logfteta <- function(eta, y) {
  n <- length(y)
  n1 <- ceiling(n / 2)
  n2 <- floor(n / 2)
  y1 <- y[1:n1]
  y2 <- y[(n1 + 1):(n1 + n2)]
  eta1 <- eta[1]
  eta2 <- eta[2]
  sum(y1) * eta1 - (length(y1) + 1) * exp(eta1) - sum(lgamma(y1+1)) + eta1 +
    sum(y2) * eta2 - (length(y2) + 1) * exp(eta2) - sum(lgamma(y2+1)) + eta2
}

n1 <- 5
n2 <- 5
n <- n1 + n2
y1 <- rpois(n1, 5)
y2 <- rpois(n2, 5)

truemode <- c(log((sum(y1) + 1) / (length(y1) + 1)), log((sum(y2) + 1) / (length(y2) + 1)))
truelogint <- function(y) lgamma(1 + sum(y)) - (1 + sum(y)) * log(length(y) + 1) - sum(lgamma(y + 1))
truelognormconst <- truelogint(y1) + truelogint(y2)

obj <- function(x) logfteta(x, c(y1, y2))

ff <- list(
  fn = obj,
  gr = function(x) numDeriv::grad(obj, x),
  he = function(x) numDeriv::hessian(obj, x)
)

grid <- expand.grid(
  theta1 = seq(0, 4, length.out = 400),
  theta2 = seq(0, 4, length.out = 400)
)

ground_truth <- data.frame(grid) %>%
  rowwise() %>%
  mutate(pdf = exp(obj(x = c(theta1, theta2))))

plot <- ggplot(ground_truth, aes(x = theta1, y = theta2, z = pdf)) +
  geom_contour(col = multi.utils::cbpalette()[1]) +
  coord_fixed(xlim = c(0.8, 2), ylim = c(0.8, 2), ratio = 1) +
  labs(x = "", y = "") +
  theme_minimal()

plot

#' Optimise using BFSG:
opt_bfgs <- aghq::optimize_theta(ff, c(1.5, 1.5), control = default_control(method = "BFGS"))
(m <- opt_bfgs$mode)
(C <- Matrix::forceSymmetric(solve(opt_bfgs$hessian)))
eigen(C)

grid_bfgs_1 <- mvQuad::createNIGrid(dim = 2, type = "GHe", level = 1)
mvQuad::rescale(grid_bfgs_1, m = m, C = C)

grid_bfgs_3 <- mvQuad::createNIGrid(dim = 2, type = "GHe", level = 3)
mvQuad::rescale(grid_bfgs_3, m = m, C = C)

grid_bfgs_5 <- mvQuad::createNIGrid(dim = 2, type = "GHe", level = 5)
mvQuad::rescale(grid_bfgs_5, m = m, C = C)

grid_bfgs_7 <- mvQuad::createNIGrid(dim = 2, type = "GHe", level = 7)
mvQuad::rescale(grid_bfgs_7, m = m, C = C)

norm_bfgs_1 <- aghq::normalize_logpost(opt_bfgs, 1, 1)
norm_bfgs_3 <- aghq::normalize_logpost(opt_bfgs, 3, 1)
norm_bfgs_5 <- aghq::normalize_logpost(opt_bfgs, 5, 1)
norm_bfgs_7 <- aghq::normalize_logpost(opt_bfgs, 7, 1)
norm_bfgs <- list(norm_bfgs_1, norm_bfgs_3, norm_bfgs_5, norm_bfgs_7)

plot_points(grid_bfgs_3)

pca_grid_bfgs_1 <- pca_rescale(m = m, C = C, s = 1, k = 1)
pca_grid_bfgs_3 <- pca_rescale(m = m, C = C, s = 1, k = 3)
pca_grid_bfgs_5 <- pca_rescale(m = m, C = C, s = 1, k = 5)
pca_grid_bfgs_7 <- pca_rescale(m = m, C = C, s = 1, k = 7)

norm_pca_bfgs_1 <- local_normalize_logpost(opt_bfgs, basegrid = pca_grid_bfgs_1, adapt = FALSE)
norm_pca_bfgs_3 <- local_normalize_logpost(opt_bfgs, basegrid = pca_grid_bfgs_3, adapt = FALSE)
norm_pca_bfgs_5 <- local_normalize_logpost(opt_bfgs, basegrid = pca_grid_bfgs_5, adapt = FALSE)
norm_pca_bfgs_7 <- local_normalize_logpost(opt_bfgs, basegrid = pca_grid_bfgs_7, adapt = FALSE)
norm_pca_bfgs <- list(norm_pca_bfgs_1, norm_pca_bfgs_3, norm_pca_bfgs_5, norm_pca_bfgs_7)

plot_points(pca_grid_bfgs_3)

pca2_grid_1 <- mvQuad::createNIGrid(dim = 2, type = "GHe", level = c(1, 1))
mvQuad::rescale(pca2_grid_1, m = m, C = C, dec.type = 1)

pca2_grid_3 <- mvQuad::createNIGrid(dim = 2, type = "GHe", level = c(3, 1))
mvQuad::rescale(pca2_grid_3, m = m, C = C, dec.type = 1)

pca2_grid_5 <- mvQuad::createNIGrid(dim = 2, type = "GHe", level = c(5, 1))
mvQuad::rescale(pca2_grid_5, m = m, C = C, dec.type = 1)

pca2_grid_7 <- mvQuad::createNIGrid(dim = 2, type = "GHe", level = c(7, 1))
mvQuad::rescale(pca2_grid_7, m = m, C = C, dec.type = 1)

norm_pca2_bfgs_1 <- local_normalize_logpost(opt_bfgs, basegrid = pca2_grid_1, adapt = FALSE)
norm_pca2_bfgs_3 <- local_normalize_logpost(opt_bfgs, basegrid = pca2_grid_3, adapt = FALSE)
norm_pca2_bfgs_5 <- local_normalize_logpost(opt_bfgs, basegrid = pca2_grid_5, adapt = FALSE)
norm_pca2_bfgs_7 <- local_normalize_logpost(opt_bfgs, basegrid = pca2_grid_7, adapt = FALSE)
norm_pca2_bfgs <- list(norm_pca2_bfgs_1, norm_pca2_bfgs_3, norm_pca2_bfgs_5, norm_pca2_bfgs_7)

plot_points(pca2_grid_3)

results <- data.frame(
  "method" = c("Truth", paste0("aghq, k = ", c(1, 3, 5, 7)), paste0("aghq-pca, k = ", c(1, 3, 5, 7)), paste0("aghq-pca2, k = ", c(1, 3, 5, 7))),
  "lognormconst" = c(truelognormconst, sapply(c(norm_bfgs, norm_pca_bfgs, norm_pca2_bfgs), function(x) x$lognormconst)),
  "type" = as.factor(c("Truth", rep("AGHQ", 4), rep("AGHQ-PCA", 4), rep("AGHQ-PCA2", 4)))
)

ggplot(results, aes(x = method, y = lognormconst, col = type)) +
  geom_point(size = 1) +
  labs(x = "Method", y = "Log normalising constant", col = "") +
  scale_color_manual(values = multi.utils::cbpalette()) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 70, vjust = 0.5, hjust = 0.5))