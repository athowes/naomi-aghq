#' Return the maximum distance from zero contained in a vector
#'
#' @param x A vector
#' @return The element with maximum absolute value
absmax <- function(x) x[which.max(abs(x))]

#' Create a facetted histogram and ECDF plots using samples from each of the four methods
#'
#' @param par The name of a parameter (latent field, hyperparameter, or model output)
#' @param i The index of the parameter, if it has dimension greater than one
#' @param return_df Should the dataframe used to create these plots be returned too?
#' @return A `ggplot2` object, or list containing a `ggplot2` object and a dataframe
histogram_and_ecdf <- function(par, i = NULL, return_df = FALSE) {
  colours <- c("#56B4E9", "#009E73", "#E69F00", "#CC79A7")

  if(!is.null(i)) {
    par_name <- paste0(par, "[", i, "]")

    df_compare <- rbind(
      data.frame(method = "TMB", samples = as.numeric(tmb$fit$sample[[par]][i, ])),
      data.frame(method = "aghq", samples = as.numeric(aghq$quad$sample[[par]][i, ])),
      data.frame(method = "tmbstan", samples = as.numeric(tmbstan$mcmc$sample[[par]][i, ]))
    )

    rhat <- signif(rstan::Rhat(tmbstan$mcmc$sample[[par]][i, ]), 3)
    ess <- signif(rstan::ess_bulk(tmbstan$mcmc$sample[[par]][i, ]), 3)

  } else {
    par_name <- paste0(par)

    df_compare <- rbind(
      data.frame(method = "TMB", samples = as.numeric(tmb$fit$sample[[par]])),
      data.frame(method = "aghq", samples = as.numeric(aghq$quad$sample[[par]])),
      data.frame(method = "tmbstan", samples = as.numeric(tmbstan$mcmc$sample[[par]]))
    )

    rhat <- signif(rstan::Rhat(tmbstan$mcmc$sample[[par]]), 3)
    ess <- signif(rstan::ess_bulk(tmbstan$mcmc$sample[[par]]), 3)
  }

  df_compare$method <- factor(df_compare$method, levels = c("TMB", "aghq", "tmbstan"))

  mean <- signif(mean(filter(df_compare, method == "tmbstan")$samples), digits = 3)
  sd <- signif(sd(filter(df_compare, method == "tmbstan")$samples), digits = 3)

  histogram_plot <- ggplot(df_compare, aes(x = samples, fill = method, col = method)) +
    geom_histogram(aes(y = after_stat(density)), alpha = 0.5, position = "identity", bins = 30) +
    theme_minimal() +
    facet_grid(method~.) +
    labs(y = "Density", fill = "Method") +
    scale_color_manual(values = colours) +
    scale_fill_manual(values = colours) +
    theme(legend.position = "none") +
    labs(title = par_name, subtitle = paste0("tmbstan: Rhat = ", rhat, ", ESS = ", ess), x = "")

  grid <- seq(from = min(df_compare$samples), to = max(df_compare$samples), length.out = 1000)

  tmb_ecdf <- stats::ecdf(filter(df_compare, method == "TMB") %>% pull(samples))
  tmb_ecdf_df <- data.frame(x = grid, ecdf = tmb_ecdf(grid), method = "TMB")

  aghq_ecdf <- stats::ecdf(filter(df_compare, method == "aghq") %>% pull(samples))
  aghq_ecdf_df <- data.frame(x = grid, ecdf = aghq_ecdf(grid), method = "aghq")

  tmbstan_ecdf <- stats::ecdf(filter(df_compare, method == "tmbstan") %>% pull(samples))
  tmbstan_ecdf_df <- data.frame(x = grid, ecdf = tmbstan_ecdf(grid), method = "tmbstan")

  # Add ECDF differences
  tmb_ecdf_df$ecdf_diff <- tmbstan_ecdf_df$ecdf - tmb_ecdf_df$ecdf
  aghq_ecdf_df$ecdf_diff <- tmbstan_ecdf_df$ecdf - aghq_ecdf_df$ecdf
  tmbstan_ecdf_df$ecdf_diff <- 0

  ks_tmb <- absmax(tmb_ecdf_df$ecdf_diff)
  ks_aghq <- absmax(aghq_ecdf_df$ecdf_diff)

  ecdf_df <- bind_rows(tmb_ecdf_df, aghq_ecdf_df, tmbstan_ecdf_df)

  ecdf_df$method <- factor(ecdf_df$method, levels = c("TMB", "aghq", "tmbstan"))

  ks_labeller <- function(x) toString(signif(abs(x), 2))

  ecdf_plot <- ggplot(ecdf_df, aes(x = x, y = ecdf_diff, col = method)) +
    geom_line() +
    geom_abline(intercept = ks_tmb, slope = 0, col = colours[1], linetype = "dashed", alpha = 0.8) +
    annotate("text", x = 1.1 * max(ecdf_df$x), y = ks_tmb, label = ks_labeller(ks_tmb), col = colours[1], alpha = 0.8) +
    geom_abline(intercept = ks_aghq, slope = 0, col = colours[2], linetype = "dashed", alpha = 0.8) +
    annotate("text", x = 1.1 * max(ecdf_df$x), y = ks_aghq, label = ks_labeller(ks_aghq), col = colours[2], alpha = 0.8) +
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

histogram_and_ecdf_list <- function(par) lapply(1:sum(names(tmb$fit$obj$env$par) == par), histogram_and_ecdf, par = par)

#' Create a dataframe of KS test statistics using samples from TMB, aghq and tmbstan
#'
#' @param par The name of a parameter (latent field, hyperparameter, or model output)
#' @param starts_with Should all parameters starting with `par` be kept? Only applicable
#' to latent field and hyperparameters
#' @param outputs Is the parameter a model output?
#' @param id Which model outputs should be kept? Usually this will be obtained from the
#' Naomi object `mf_out` and used to only select model outputs at the finest level
#' @return A dataframe
to_ks_df <- function(par, starts_with = FALSE, outputs = FALSE, id = NULL) {
  par_names <- par
  all_par_names <- names(tmb$fit$obj$env$par)
  if(starts_with) par_names <- all_par_names[stringr::str_starts(all_par_names, par)]
  unique_par_names <- unique(par_names)

  # If outputs = TRUE then we use id to subset the samples to only to outputs at
  # the finest level of granularity, avoiding double counting aggregate measures
  process_samples <- function(x, outputs = FALSE, id = NULL) {
    samples <- x[unique_par_names]
    samples <- lapply(samples, function(x) {
      if(outputs) x <- x[id, ]
      as.data.frame(t(x))
    })
    if(!outputs) {
      table <- table(all_par_names)
      for(par in unique_par_names) {
        par_length <- table[par]
        if(par_length > 1) {
          par_colnames <- paste0(par, "[", 1:par_length, "]")
        } else {
          par_colnames <- paste0(par)
        }
        colnames(samples[[par]]) <- par_colnames
      }
    }
    dplyr::bind_cols(samples)
  }

  samples_tmb <- process_samples(tmb$fit$sample, outputs = outputs, id = id)
  samples_aghq <- process_samples(aghq$quad$sample, outputs = outputs, id = id)
  samples_tmbstan <- process_samples(tmbstan$mcmc$sample, outputs = outputs, id = id)

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

#' Create a scatterplot and jitterplot of the KS test statistics
#'
#' @param ks_df The output of `to_ks_df`
#' @param par Parameter name (only used for labelling)
#' @param method1 Samples from this method will be used as the first entry in the KS test
#' @param method2 Samples from this method will be used as the second entry in the KS test
#' @return A `ggplot2` object
ks_plot <- function(ks_df, par, method1 = "TMB", method2 = "aghq", alpha = 0.5) {
  wide_ks_df <- pivot_wider(ks_df, names_from = "method", values_from = "ks")
  wide_ks_df$ks_diff <- wide_ks_df[[method1]] - wide_ks_df[[method2]]
  mean_ks_diff <- mean(wide_ks_df$ks_diff)

  jitterplot <- wide_ks_df %>%
    ggplot(aes(x = "", y = ks_diff)) +
    geom_jitter(shape = 1, width = 0.05, height = 0, alpha = alpha) +
    scale_y_continuous(labels = scales::percent) +
    geom_boxplot(fill = NA, width = 0.2, outlier.shape = NA) +
    labs(
      caption = paste0(
        "KS tests for ", par, " of length ", nrow(ks_df) / length(unique(ks_df$method)),
        " with a mean KS difference of ", 100 * signif(mean_ks_diff, 3), "%.\n",
        "(If >0% then ", method1, " more different to tmbstan, and if <0% then ", method2, " more different)"),
      x = "", y = paste0("KS(", method1, ", tmbstan) - KS(", method2,", tmbstan)")
    ) +
    theme_minimal()

  scatterplot <- ggplot() +
    geom_point(data = wide_ks_df, aes(x = .data[[method1]], y = .data[[method2]]), shape = 1, alpha = alpha) +
    scale_x_continuous(labels = scales::percent) +
    scale_y_continuous(labels = scales::percent) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    labs(
      x = paste0("KS(", method1, ", tmbstan)"), y = paste0("KS(", method2,", tmbstan)"),
      fill = "KS difference"
    ) +
    theme_minimal() +
    guides(fill = "none")

  scatterplot + jitterplot +
    plot_layout(widths = c(2, 1))
}

#' Create a density plot and ridgeplot of the KS test statistics
#'
#' @param ks_df The output of `to_ks_df`
#' @param par Parameter name (only used for labelling)
#' @param method1 Samples from this method will be used as the first entry in the KS test
#' @param method2 Samples from this method will be used as the second entry in the KS test
#' @return A `ggplot2` object
ks_plot_many <- function(ks_summary, method1, method2) {
  ks_method1 <- paste0("KS(", method1, ", tmbstan)")
  ks_method2 <- paste0("KS(", method2, ", tmbstan)")

  xy_length <- min(1, max(ks_summary[[ks_method1]], ks_summary[[ks_method2]]) + 0.03)

  densityplot <- ggplot() +
    stat_density_2d(data = ks_summary, aes(x = .data[[ks_method1]], y = .data[[ks_method2]], linetype = type), col = "black") +
    xlim(0, xy_length) +
    ylim(0, xy_length) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    labs(x = ks_method1, y = ks_method2) +
    theme_minimal() +
    guides(fill = "none") +
    theme(legend.position = "bottom")

  ks_summary[["KS difference"]] <- ks_summary[[ks_method1]] - ks_summary[[ks_method2]]

  ridgeplot <- ggplot(ks_summary, aes(y = type, x = `KS difference`)) +
    ggridges::geom_density_ridges(alpha = 0.7, fill = NA, aes(linetype = type)) +
    coord_flip() +
    scale_linetype_manual(values = c("solid", "dashed")) +
    labs(y = "", x = paste0(ks_method1, " - ", ks_method2)) +
    guides(linetype = "none") +
    theme_minimal()

  densityplot + ridgeplot
}
