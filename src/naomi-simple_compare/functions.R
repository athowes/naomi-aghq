absmax <- function(x) x[which.max(abs(x))]

#' Create a facetted histogram plot of the named parameter using samples from each of the four methods
#' As well as an ECDF plot
histogram_and_ecdf <- function(par, i = NULL, return_df = FALSE) {
  colours <- c("#56B4E9", "#009E73", "#E69F00", "#CC79A7")

  if(!is.null(i)) {
    par_name <- paste0(par, "[", i, "]")

    df_compare <- rbind(
      data.frame(method = "TMB", samples = as.numeric(tmb$fit$sample[[par]][i, ])),
      data.frame(method = "aghq", samples = as.numeric(aghq$quad$sample[[par]][i, ])),
      data.frame(method = "adam", samples = as.numeric(adam$sample[[par]][i, ])),
      data.frame(method = "tmbstan", samples = as.numeric(unlist(rstan::extract(tmbstan$mcmc, pars = par)[[par]][, i])))
    )
  } else {
    par_name <- paste0(par)

    df_compare <- rbind(
      data.frame(method = "TMB", samples = as.numeric(tmb$fit$sample[[par]])),
      data.frame(method = "aghq", samples = as.numeric(aghq$quad$sample[[par]])),
      data.frame(method = "adam", samples = as.numeric(adam$sample[[par]])),
      data.frame(method = "tmbstan", samples = as.numeric(unlist(rstan::extract(tmbstan$mcmc, pars = par))))
    )
  }

  df_compare$method <- factor(df_compare$method, levels = c("TMB", "aghq", "adam", "tmbstan"))

  mean <- df_compare %>%
    filter(method == "tmbstan") %>%
    summarise(mean = mean(samples)) %>%
    pull(mean) %>%
    round(digits = 3)

  sd <- df_compare %>%
    filter(method == "tmbstan") %>%
    summarise(sd = sd(samples)) %>%
    pull(sd) %>%
    round(digits = 3)

  histogram_plot <- ggplot(df_compare, aes(x = samples, fill = method, col = method)) +
    geom_histogram(aes(y = after_stat(density)), alpha = 0.5, position = "identity", bins = 30) +
    theme_minimal() +
    facet_grid(method~.) +
    labs(y = "Density", fill = "Method") +
    scale_color_manual(values = colours) +
    scale_fill_manual(values = colours) +
    theme(legend.position = "none") +
    labs(title = par_name, subtitle = paste0("Mean = ", mean, ", SD = ", sd, " (from tmbstan)"))

  grid <- seq(from = min(df_compare$samples), to = max(df_compare$samples), length.out = 1000)

  tmb_ecdf <- stats::ecdf(filter(df_compare, method == "TMB") %>% pull(samples))
  tmb_ecdf_df <- data.frame(x = grid, ecdf = tmb_ecdf(grid), method = "TMB")

  aghq_ecdf <- stats::ecdf(filter(df_compare, method == "aghq") %>% pull(samples))
  aghq_ecdf_df <- data.frame(x = grid, ecdf = aghq_ecdf(grid), method = "aghq")

  adam_ecdf <- stats::ecdf(filter(df_compare, method == "adam") %>% pull(samples))
  adam_ecdf_df <- data.frame(x = grid, ecdf = adam_ecdf(grid), method = "adam")

  tmbstan_ecdf <- stats::ecdf(filter(df_compare, method == "tmbstan") %>% pull(samples))
  tmbstan_ecdf_df <- data.frame(x = grid, ecdf = tmbstan_ecdf(grid), method = "tmbstan")

  # Add ECDF differences
  tmb_ecdf_df$ecdf_diff <- tmbstan_ecdf_df$ecdf - tmb_ecdf_df$ecdf
  aghq_ecdf_df$ecdf_diff <- tmbstan_ecdf_df$ecdf - aghq_ecdf_df$ecdf
  adam_ecdf_df$ecdf_diff <- tmbstan_ecdf_df$ecdf - adam_ecdf_df$ecdf
  tmbstan_ecdf_df$ecdf_diff <- 0

  ks_tmb <- absmax(tmb_ecdf_df$ecdf_diff)
  ks_aghq <- absmax(aghq_ecdf_df$ecdf_diff)
  ks_adam <- absmax(adam_ecdf_df$ecdf_diff)

  ecdf_df <- bind_rows(tmb_ecdf_df, aghq_ecdf_df, adam_ecdf_df, tmbstan_ecdf_df)

  ecdf_df$method <- factor(ecdf_df$method, levels = c("TMB", "aghq", "adam", "tmbstan"))

  ks_labeller <- function(x) toString(round(abs(x), 2))

  ecdf_plot <- ggplot(ecdf_df, aes(x = x, y = ecdf_diff, col = method)) +
    geom_line() +
    geom_abline(intercept = ks_tmb, slope = 0, col = colours[1], linetype = "dashed", alpha = 0.8) +
    annotate("text", x = 1.1 * max(ecdf_df$x), y = ks_tmb, label = ks_labeller(ks_tmb), col = colours[1], alpha = 0.8) +
    geom_abline(intercept = ks_aghq, slope = 0, col = colours[2], linetype = "dashed", alpha = 0.8) +
    annotate("text", x = 1.1 * max(ecdf_df$x), y = ks_aghq, label = ks_labeller(ks_aghq), col = colours[2], alpha = 0.8) +
    geom_abline(intercept = ks_adam, slope = 0, col = colours[3], linetype = "dashed", alpha = 0.8) +
    annotate("text", x = 1.1 * max(ecdf_df$x), y = ks_adam, label = ks_labeller(ks_adam), col = colours[3], alpha = 0.8) +
    scale_color_manual(values = colours) +
    labs(x = "", y = "ECDF difference") +
    guides(col = "none") +
    coord_cartesian(xlim = c(min(ecdf_df$x), max(ecdf_df$x)), clip = "off") +
    theme_minimal() +
    theme(plot.margin = unit(c(1, 3, 1, 1), "lines"))

  plot <- histogram_plot + ecdf_plot

  if(return_df) {
    return(list(plot = plot, df = df_compare))
  } else {
    plot
  }
}

#' Create a dataframe of samples from TMB, aghq, adam and tmbstan for any parameters starting with par
to_ks_df <- function(par, starts_with = FALSE) {
  par_names <- par

  if(starts_with) {
    all_par_names <- names(tmb$fit$obj$env$par)
    par_names <- all_par_names[stringr::str_starts(all_par_names, par)]
  }

  unique_par_names <- unique(par_names)

  samples_tmb <- tmb$fit$sample[unique_par_names]
  samples_tmb <- lapply(samples_tmb, function(x) as.data.frame(t(x)))

  samples_aghq <- aghq$quad$sample[unique_par_names]
  samples_aghq <- lapply(samples_aghq, function(x) as.data.frame(t(x)))

  samples_adam <- adam$sample[unique_par_names]
  samples_adam <- lapply(samples_adam, function(x) as.data.frame(t(x)))

  samples_tmbstan <- rstan::extract(tmbstan$mcmc, pars = unique_par_names)
  samples_tmbstan <- lapply(samples_tmbstan, function(x) as.data.frame(x))

  table <- table(par_names)
  for(par in unique_par_names) {
    par_length <- table[par]
    if(par_length > 1) {
      par_colnames <- paste0(par, "[", 1:par_length, "]")
    } else {
      par_colnames <- paste0(par)
    }
    colnames(samples_tmb[[par]]) <- par_colnames
    colnames(samples_aghq[[par]]) <- par_colnames
    colnames(samples_adam[[par]]) <- par_colnames
    colnames(samples_tmbstan[[par]]) <- par_colnames
  }

  samples_tmb <- dplyr::bind_cols(samples_tmb)
  samples_aghq <- dplyr::bind_cols(samples_aghq)
  samples_adam <- dplyr::bind_cols(samples_adam)
  samples_tmbstan <- dplyr::bind_cols(samples_tmbstan)

  n <- ncol(samples_tmbstan)
  ks_tmb <- numeric(n)
  ks_aghq <- numeric(n)
  ks_adam <- numeric(n)
  for(i in 1:n) {
    ks_tmb[i] <- inf.utils::ks_test(samples_tmb[, i], samples_tmbstan[, i])$D
    ks_aghq[i] <- inf.utils::ks_test(samples_aghq[, i], samples_tmbstan[, i])$D
    ks_adam[i] <- inf.utils::ks_test(samples_adam[, i], samples_tmbstan[, i])$D
  }

  rbind(
    data.frame(method = "TMB", ks = ks_tmb, par = names(samples_tmbstan), index = 1:n),
    data.frame(method = "aghq", ks = ks_aghq, par = names(samples_tmbstan), index = 1:n),
    data.frame(method = "adam", ks = ks_adam, par = names(samples_tmbstan), index = 1:n)
  )
}

ks_plot <- function(ks_df, par, method1 = "TMB", method2 = "aghq") {
  wide_ks_df <- pivot_wider(ks_df, names_from = "method", values_from = "ks")
  wide_ks_df$ks_diff <- wide_ks_df[[method1]] - wide_ks_df[[method2]]
  mean_ks_diff <- mean(wide_ks_df$ks_diff)

  jitterplot <- wide_ks_df %>%
    ggplot(aes(x = "", y = ks_diff)) +
    geom_jitter(width = 0.1, alpha = 0.5) +
    labs(
      title = paste0("Mean KS difference* is ", round(mean_ks_diff, 3)),
      caption = paste0("*If >0 then ", method1, " more different to tmbstan, and if <0 then ", method2, " more different"),
      x = "", y = paste0("KS(", method1, ", tmbstan) - KS(", method2,", tmbstan)")
    ) +
    theme_minimal()

  xy_length <- min(1, max(wide_ks_df[[method1]], wide_ks_df[[method2]]) + 0.05)

  scatterplot <- ggplot(wide_ks_df, aes(x = .data[[method1]], y = .data[[method2]])) +
    geom_point(alpha = 0.5) +
    annotate(geom = "polygon", x = c(-Inf, Inf, Inf), y = c(-Inf, Inf, -Inf), fill = multi.utils::cbpalette()[2], alpha = 0.2) +
    annotate(geom = "polygon", x = c(-Inf, Inf, -Inf), y = c(-Inf, Inf, Inf), fill = multi.utils::cbpalette()[6], alpha = 0.2) +
    xlim(0, xy_length) +
    ylim(0, xy_length) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    labs(
      title = paste0("KS tests for ", par, " of length ", max(ks_df$index)),
      x = paste0("KS(", method1, ", tmbstan)"), y = paste0("KS(", method2,", tmbstan)")
    ) +
    theme_minimal()

  scatterplot + jitterplot
}
