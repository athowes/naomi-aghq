library(patchwork)
library(Matrix)
library(data.table)

#' Test code for inverting Hessian here

quad <- readRDS("depends/aghq.rds")$quad

#' Fit with k = 1 so there is only one modesandhessians entry
nrow(quad$modesandhessians)

H <- quad$modesandhessians$H[[1]]

H_df <- reshape2::melt(as.matrix(H))

pdf("H-matrix.pdf", h = 7, w = 6.25)

histogram <- ggplot(H_df, aes(x = log10(value))) +
  geom_histogram(alpha = 0.8) +
  labs(x = "log10(H[i, j])", y = "Count") +
  theme_minimal()

heatmap <- H_df %>%
  mutate(value_na = ifelse(value == 0, NA, value)) %>%
  ggplot(aes(x = Var1, y = Var2, fill = log10(value))) +
  labs(x = "", y = "", fill = "log10(H[i, j])") +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(x = "i", y = "j") +
  theme_minimal()

heatmap / histogram

dev.off()

saveRDS(H, "H.rds")

#' Aiming to get the marginal standard deviations without inverting the whole Hessian
diag_of_inv_naive <- function(H) {
  diag(solve(H))
}

# Stolen from Patrick Brown
diag_of_inv_patrick <- function(H) {
  L <- Matrix::expand(Matrix::Cholesky(H, LDL = FALSE, super = FALSE))
  L$Linv <- Matrix::solve(L$L)
  L$LinvDf <- data.table::data.table(
    col = rep(1:nrow(L$Linv), diff(L$Linv@p)),
    x = L$Linv@x^2
  )
  varDiag <- L$LinvDf[, .(sum = sum(x)), by = col]
  varDiagMat <- Diagonal(nrow(varDiag), varDiag$sum)
  varDiagMatP <- crossprod(L$P, varDiagMat) %*% L$P
  varDiagMatP@x
}

microbenchmark::microbenchmark(
  diag_of_inv_naive(H),
  diag_of_inv_patrick(H)
)

ggplot() +
  geom_point(aes(x = diag_of_inv_naive(H), y = diag_of_inv_patrick(H))) +
  labs(x = "Naive", y = "Patrick") +
  theme_minimal()

stopifnot(abs(diag_of_inv_naive(H) - diag_of_inv_patrick(H)) < 10E-12)
