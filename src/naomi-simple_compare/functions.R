#' Create a facetted histogram plot of the named parameter using samples from each of the three methods
histogram_plot <- function(par) {
  df_compare <- rbind(
    data.frame(method = "TMB", samples = as.numeric(tmb$fit$sample[[par]])),
    data.frame(method = "aghq", samples = as.numeric(aghq$quad$sample[[par]])),
    data.frame(method = "tmbstan", samples = as.numeric(unlist(rstan::extract(tmbstan$mcmc, pars = par))))
  )

  df_compare %>%
    group_by(method) %>%
    summarise(n = n())

  ggplot(df_compare, aes(x = samples, fill = method, col = method)) +
    geom_histogram(aes(y = after_stat(density)), alpha = 0.5, position = "identity", bins = 30) +
    theme_minimal() +
    facet_grid(method~.) +
    labs(x = paste0(par), y = "Density", fill = "Method") +
    scale_color_manual(values = multi.utils::cbpalette()) +
    scale_fill_manual(values = multi.utils::cbpalette()) +
    theme(legend.position = "none") +
    labs(title = paste0(par))
}

#' Create a dataframe of samples from TMB, aghq and tmbstan for any parameters starting with par
to_ks_df <- function(par) {
  all_par_names <- names(tmb$fit$obj$env$par)
  par_names <- all_par_names[stringr::str_starts(all_par_names, par)]
  unique_par_names <- unique(par_names)

  samples_tmb <- tmb$fit$sample[unique_par_names]
  samples_tmb <- lapply(samples_tmb, function(x) as.data.frame(t(x)))

  samples_aghq <- aghq$quad$sample[unique_par_names]
  samples_aghq <- lapply(samples_aghq, function(x) as.data.frame(t(x)))

  samples_tmbstan <- rstan::extract(tmbstan$mcmc, pars = unique_par_names)
  samples_tmbstan <- lapply(samples_tmbstan, function(x) as.data.frame(x))

  table <- table(par_names)
  unique_par_names <- unique(par_names)
  for(par in unique_par_names) {
    par_length <- table[par]
    if(par_length > 1) {
      par_colnames <- paste0(par, "[", 1:par_length, "]")
    } else {
      par_colnames <- paste0(par)
    }
    colnames(samples_tmb[[par]]) <- par_colnames
    colnames(samples_aghq[[par]]) <- par_colnames
    colnames(samples_tmbstan[[par]]) <- par_colnames
  }

  samples_tmb <- dplyr::bind_cols(samples_tmb)
  samples_aghq <- dplyr::bind_cols(samples_aghq)
  samples_tmbstan <- dplyr::bind_cols(samples_tmbstan)

  n <- ncol(samples_tmbstan)
  ks_tmb <- numeric(n)
  ks_aghq <- numeric(n)
  for(i in 1:n) {
    ks_tmb[i] <- inf.utils::ks_test(samples_tmb[, i], samples_tmbstan[, i])$D
    ks_aghq[i] <- inf.utils::ks_test(samples_aghq[, i], samples_tmbstan[, i])$D
  }

  rbind(
    data.frame(method = "TMB", ks = ks_tmb, par = names(samples_tmbstan), index = 1:n),
    data.frame(method = "aghq", ks = ks_aghq, par = names(samples_tmbstan), index = 1:n)
  )
}


ks_plot <- function(ks_df, par) {
  wide_ks_df <- pivot_wider(ks_df, names_from = "method", values_from = "ks") %>%
    mutate(ks_diff = TMB - aghq)

  mean_ks_diff <- mean(wide_ks_df$ks_diff)

  boxplot <- wide_ks_df %>%
    ggplot(aes(x = ks_diff)) +
    geom_boxplot(width = 0.5) +
    coord_flip() +
    labs(
      title = paste0("Mean KS difference is ", mean_ks_diff),
      subtitle = ">0 then TMB more different to tmbstan, <0 then aghq more different",
      x = "KS(TMB, tmbstan) - KS(aghq, tmbstan)"
    ) +
    theme_minimal()

  xy_length <- min(1, max(wide_ks_df$TMB, wide_ks_df$aghq) + 0.05)

  scatterplot <- ggplot(wide_ks_df, aes(x = TMB, y = aghq)) +
    geom_point(alpha = 0.5) +
    xlim(0, xy_length) +
    ylim(0, xy_length) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    labs(
      title = paste0("KS tests for ", par, " of length ", max(ks_df$index)),
      subtitle = "Values along y = x have similar KS",
      x = "KS(aghq, tmbstan)", y = "KS(TMB, tmbstan)"
    ) +
    theme_minimal()

  cowplot::plot_grid(scatterplot, boxplot, ncol = 2, rel_widths = c(1.3, 1))
}
