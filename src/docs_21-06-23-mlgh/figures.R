TMB::compile("2d.cpp")
dyn.load(TMB::dynlib("2d"))

obj <- TMB::MakeADFun(data = list(), parameters = list(theta1 = 0, theta2 = 0), DLL = "2d")

grid <- expand.grid(
  theta1 = seq(-8, 8, length.out = 400),
  theta2 = seq(-8, 8, length.out = 400)
)

ground_truth <- cbind(grid, pdf = apply(grid, 1, function(x) exp(-1 * obj$fn(x))))

opt <- nlminb(
  start = obj$par,
  objective = obj$fn,
  gradient = obj$gr,
  control = list(iter.max = 1000, trace = 0)
)

sd_out <- TMB::sdreport(
  obj,
  par.fixed = opt$par,
  getJointPrecision = TRUE
)

mu <- opt$par
cov <- sd_out$cov.fixed

figA0 <- ggplot(ground_truth, aes(x = theta1, y = theta2, z = pdf)) +
  geom_contour(col = "lightgrey") +
  coord_fixed(xlim = c(-8, 8), ylim = c(-8, 8), ratio = 1) +
  labs(x = "", y = "") +
  theme_minimal() +
  guides(size = "none") +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank()
  )

gg <- mvQuad::createNIGrid(2, "GHe", 3)

add_points <- function(figA0, gg) {

  points <- mvQuad::getNodes(gg) %>%
    as.data.frame() %>%
    mutate(weights = mvQuad::getWeights(gg))

  colnames(points) <- c("theta1", "theta2", "weights")

  figA0 +
    geom_point(
      data = points,
      aes(x = theta1, y = theta2, size = weights),
      alpha = 0.8,
      col = "#009E73",
      inherit.aes = FALSE
    ) +
    scale_size_continuous(range = c(1, 2))
}

figA1 <- add_points(figA0, gg) +
  labs(size = "")

#' Adapt by the mean
gg2 <- gg
mvQuad::rescale(gg2, m = mu, C = diag(c(1, 1)), dec.type = 1)

figA2 <- add_points(figA0, gg2) +
  labs(size = "")

#' Adapt by the lower Cholesky
gg3 <- gg
mvQuad::rescale(gg3, m = mu, C = cov, dec.type = 2)

figA3 <- add_points(figA0, gg3) +
  labs(size = "")

#' Adapt by the spectral
gg4 <- gg
mvQuad::rescale(gg4, m = mu, C = cov, dec.type = 1)

figA4 <- add_points(figA0, gg4) +
  labs(size = "")

#' PCA-AGHQ
gg5 <- mvQuad::createNIGrid(2, "GHe", level = c(3, 1))
mvQuad::rescale(gg5, m = mu, C = cov, dec.type = 1)

figA5 <- add_points(figA0, gg5) +
  labs(size = "")
