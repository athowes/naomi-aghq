#' This development script is to explore the idea of using a different value of k
#' for each hyperparameter to build a product grid which in the end is not prohibitively
#' expensive. The first method I thought of to select the value of k for a particular
#' dimension was based on the standard deviation of that hyperparameter when fitted
#' with a cheaper model e.g. TMB empirical Bayes plus Gaussian approximation.

#' Given TMB fit chose k for each dimension
hypers <- names(fit$par)

#' Should be 24
length(hypers)

means <- vector(mode = "numeric", length = length(hypers))
sds <- vector(mode = "numeric", length = length(hypers))

for(i in seq_along(hypers)) {
  means[i] <- mean(fit$sample[[hypers[i]]])
  sds[i] <- sd(fit$sample[[hypers[i]]])
}

#' Interval plot with parameter names along one side

pdf("hyper-intervals.pdf", h = 5, w = 6.25)

data.frame(name = hypers, mean = means, sd = sds) %>%
  ggplot(aes(x = name, y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean - 2 * sd, ymax = mean + 2 * sd)) +
  coord_flip() +
  labs(x = "", y = "") +
  scale_x_discrete(position = "top") +
  theme_minimal()

dev.off()

#' Total number of points if it were two grid points per dimension
formatC(2^length(hypers), format = "e")

#' Total number of points if it were three grid points per dimension
formatC(3^length(hypers), format = "e")

#' What about an alternative strategy
cut_offs <- c(0, 1, 1.5)
k_options <- c(1, 2, 3)

ks <- vector(mode = "numeric", length = length(hypers))

ks[sds > cut_offs[1]] <- k_options[1]
ks[sds > cut_offs[2]] <- k_options[2]
ks[sds > cut_offs[3]] <- k_options[3]

#' The total number of grid points in a product grid would be
prod(ks)

grids <- list()

for(i in seq_along(hypers)) {
  grids[[i]] <- mvQuad::createNIGrid(1, "GHe", level = ks[i], "product")
}

#' Is there a way to take a product of mvQuad grids outside of the function?
grids[[1]]

#' Actually: looks like there is a way to do this with mvQuad::createNIGrid
#' See the argument `level.trans`
#'
#' Alternatively `level.trans` can be a function, which takes (n x d)-matrix and returns
#' a matrix with the same dimensions (see the example; this feature is particularly useful
#' for the 'sparse' construction rule, to account for different importance of the dimensions).
#'
#' Start of examples

#' 1D grid: closed Newton-Cotes formula of degree 1 (trapeziodal-rule)
grid <- mvQuad::createNIGrid(dim = 1, type = "cNC1", level = 10)
print(grid)
plot(grid)

#' 2D grid: nested Gauss-Legendre rule
grid <- mvQuad::createNIGrid(dim = 2, type = c("GLe", "nLe"), level = c(4, 7))
mvQuad::rescale(grid, domain = rbind(c(-1, 1), c(-1, 1)))
plot(grid)
print(grid)

f <- function(x) {
  1 - x[, 1]^2 * x[, 2]^2
}

mvQuad::quadrature(f = f, grid = grid)

#' Level transformation
level_tranformation <- function(x){
  tmp <- as.matrix(x)
  tmp[, 2] <- 2 * tmp[ ,2]
  return(tmp)
}

nw <- mvQuad::createNIGrid(dim = 2, type = "cNC1", level = 3, level.trans = level_tranformation, ndConstruction = "sparse")
plot(nw)

#' Now let's try to use the `level.trans` argument ourselves
#' Note that `level` is the accuracy level a.k.a. number of grid points
#' Ah! I don't think we need this `level.trans`, we can just use the following:
#' "The argument type and level can also be vector-value, different for each dimension
#' (the later only for "product rule"; see examples)"

grid <- mvQuad::createNIGrid(dim = 2, type = "GHe", level = c(3, 5))
plot(grid)

#' Looks fine to me!
grid <- mvQuad::createNIGrid(dim = length(hypers), "GHe", level = ks)
grid

#' Create a function to make this type of grid
sd_levels_ghe_grid <- function(dim, level, cut_off, mean, sd) {
  stopifnot(length(level) == length(cut_off))
  stopifnot(dim == length(sd))
  levels <- vector(mode = "numeric", length = dim)
  for(i in seq_along(cut_off)) levels[sd > cut_offs[i]] <- level[i]
  grid <- mvQuad::createNIGrid(dim = dim, "GHe", level = levels)
  grid
}

sd_levels_ghe_grid(dim = length(hypers), level = c(1, 2), cut_off = c(0, 1.5), sd = sds)
