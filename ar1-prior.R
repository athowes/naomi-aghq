#' This is the prior being used in Naomi at the moment
nsim <- 10000
logit_phi <- rnorm(nsim, 0, sd = 2.582)
phi <- (2 * plogis(logit_phi)) - 1

data.frame(
  x = phi
) %>%
  ggplot(aes(x = x)) +
    geom_histogram() +
    labs(x = "phi", y = "p(phi)") +
    theme_minimal()

#' This is the R-INLA default
theta <- rnorm(nsim, 0, sd = 0.15)
from.theta <- function(x) 2 * exp(x) / (1 + exp(x)) - 1
phi2 <- from.theta(theta)

data.frame(
  x = phi2
) %>%
  ggplot(aes(x = x)) +
  geom_histogram() +
  labs(x = "phi", y = "p(phi)") +
  theme_minimal()
