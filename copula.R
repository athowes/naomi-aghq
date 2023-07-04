#' Y ~ N(mu, cov)
#' Y = mu + LZ
#' L is s.th LL^T = cov

rho <- 0.99
cor <- matrix(c(1, rho, rho, 1), ncol = 2, nrow = 2)
chol(cor)
