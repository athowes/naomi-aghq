simulate_multi_risk <- function(prob_Y015_019 = c(0.6, 0.25, 0.12, 0.03),
                                prob_Y020_024 = c(0.2, 0.4, 0.35, 0.05),
                                prob_Y025_029 = c(0.05, 0.5, 0.4, 0.05),
                                n = 30,
                                m = 50) {

  #' Using three age groups here
  age_groups <- c("Y015_019", "Y020_024", "Y025_029")

  #' Assume these four categories to be exhaustive and mutually exclusive
  indicators <- c("nosex12m", "sexcohab", "sexnonreg", "sexpaid12m")
  K <- length(indicators)

  df <- crossing(
    age_group = age_groups,
    area_id = 1:n,
    indicator = indicators,
    sample_id = 1:m
  )

  #' N is the number of areas
  samp_Y015_019 <- tall_rmultinomial(N = m * n, prob = prob_Y015_019)
  samp_Y020_024 <- tall_rmultinomial(N = m * n, prob = prob_Y020_024)
  samp_Y025_029 <- tall_rmultinomial(N = m * n, prob = prob_Y025_029)

  df <- df %>%
    arrange(age_group) %>%
    mutate(
      y = c(samp_Y015_019, samp_Y020_024, samp_Y025_029)
    )

  return(df)
}

#' Tall format: (1 x N x 4) length vector
tall_rmultinomial <- function(N, prob) {
  replicate(N, rmultinom(1, 1, prob = prob), simplify = TRUE) %>%
    as.vector()
}
