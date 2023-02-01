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
data.frame(name = hypers, mean = means, sd = sds) %>%
  ggplot(aes(x = name, y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean - 2 * sd, ymax = mean + 2 * sd)) +
  coord_flip() +
  labs(x = "", y = "") +
  scale_x_discrete(position = "top") +
  theme_minimal()

#' Total number of points if it were two grid points per dimension
formatC(2^length(hypers), format = "e")

#' Total number of points if it were three grid points per dimension
formatC(3^length(hypers), format = "e")

cut_offs <- c(0, 1, 1.5)
k_options <- c(1, 2, 3)

ks <- vector(mode = "numeric", length = length(hypers))

ks[sds > cut_offs[1]] <- k_options[1]
ks[sds > cut_offs[2]] <- k_options[2]
ks[sds > cut_offs[3]] <- k_options[3]

#' The total number of grid points in a product grid would be
prod(k)

grids <- list()

for(i in seq_along(hypers)) {
  grids[[i]] <- mvQuad::createNIGrid(1, "GHe", level = ks[i], "product")
}

#' Is there a way to take a product of mvQuad grids outside of the function?
grids[[1]]

#' Actually it looks like there is a way to do this natrually with mvQuad::createNIGrid
#' See the argument `level.trans`
