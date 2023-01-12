quad <- sparse_quad
M <- 10
transformation <- quad$transformation
interpolation <- "auto"

numcores <- getOption("mc.cores", 1L)

if (.Platform$OS.type == "windows") numcores <- 1

d <- dim(quad$modesandhessians$H[[1]])[1]

simlist <- quad$modesandhessians

if (numcores > 1) {
  simlist$L <- parallel::mclapply(simlist$H, function(h) Matrix::Cholesky(as(Matrix::forceSymmetric(h), "sparseMatrix"), perm = TRUE, LDL = FALSE), mc.cores = numcores)
} else {
  simlist$L <- lapply(simlist$H, function(h) Matrix::Cholesky(as(Matrix::forceSymmetric(h), "sparseMatrix"), perm = TRUE, LDL = FALSE))
}

simlist$lambda <- exp(quad$normalized_posterior$nodesandweights$logpost_normalized) * quad$normalized_posterior$nodesandweights$weights

if (M == 1) {
  k <- which(stats::rmultinom(M, 1, simlist$lambda) == 1)
} else {
  k <- apply(stats::rmultinom(M, 1, simlist$lambda), 2, function(x) which(x == 1)) #' This is where the error is
}

#' There is one particular value which is very negative
plot(simlist$lambda)

#' Which stems from the negative weight here
plot(quad$normalized_posterior$nodesandweights$weights)

#' Which can just be reproduced by
mvQuad::getWeights(mvQuad::createNIGrid(24, "GHe", 2, "sparse"))

Z <- lapply(split(matrix(stats::rnorm(M * d), nrow = M), k), matrix, nrow = d)

samps <- mapply(
  function(.x, .y) as.numeric(Matrix::solve(simlist$L[[as.numeric(.y)]], Matrix::solve(simlist$L[[as.numeric(.y)]], .x, system = "Lt"), system = "Pt")) + do.call(cbind, rep(list(simlist$mode[[as.numeric(.y)]]), ncol(.x))),
  Z,
  names(Z),
  SIMPLIFY = FALSE
)

ord <- numeric(length(k))
cumtab <- cumsum(c(0, table(k)))
cumtab <- cumtab[-length(cumtab)]
cnt <- numeric(length(unique(k)))
names(cnt) <- sort(unique(k))

for (i in 1:length(k)) {
  wc <- which(names(cnt) == k[i])
  cnt[wc] <- cnt[wc] + 1
  ord[i] <- cumtab[wc] + cnt[wc]
}

samps <- Reduce(cbind, samps)
samps <- samps[, ord]
md <- length(quad$optresults$mode)
thetanames <- colnames(quad$normalized_posterior$nodesandweights)[1:md]
theta <- simlist[k, thetanames]

if (!is.matrix(samps)) {
  samps <- rbind(samps)
  rownames(samps) <- NULL
}

if (!inherits(theta, "data.frame")) theta <- data.frame(theta1 = theta)

out <- list(samps = samps, theta = theta)
class(quad) <- "aghq"
out$thetasamples <- sample_marginal(quad, M, transformation, interpolation)

out
