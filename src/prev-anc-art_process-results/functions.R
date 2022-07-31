draw_boxplots <- function(results) {
  map(results, "comparison_results") %>%
    bind_rows(.id = "sim_id") %>%
    pivot_longer(cols = c("mean", "sd"), names_to = "type", values_to = "value") %>%
    left_join(true_values, by = "parameter") %>%
    mutate(true_value = ifelse(type == "SD", NA, true_value)) %>%
    ggplot(aes(x = method, y = value, fill = method)) +
      geom_boxplot() +
      geom_hline(aes(yintercept = true_value), linetype = "dashed") +
      facet_wrap(parameter ~ type, scales = "free", ncol = 2) +
      scale_fill_manual(values = cbpalette) +
      labs(x = "", y = "Estimate", fill = "Inference method") +
      theme_minimal()
}

draw_scatterplots <- function(results) {
  df <- map(results, "comparison_results") %>%
    bind_rows(.id = "sim_id") %>%
    pivot_longer(cols = c("mean", "sd"), names_to = "type", values_to = "value") %>%
    left_join(true_values, by = "parameter") %>%
    mutate(true_value = ifelse(type == "SD", NA, true_value))

  df <- df %>%
    filter(method != "tmbstan") %>%
    left_join(
      df %>%
        filter(method == "tmbstan") %>%
        rename(tmbstan_value = value) %>%
        select(sim_id, type, parameter, tmbstan_value),
      by = c("sim_id", "type", "parameter"),
    )

  ggplot(df, aes(x = tmbstan_value, y = value, col = method)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    facet_wrap(type ~ parameter) +
    scale_color_manual(values = cbpalette) +
    labs(x = "Gold standard (tmbstan)", y = "Value", col = "Inference method") +
    theme_minimal() +
    theme(
      legend.position = "bottom"
    )
}

draw_ksplots_D <- function(results) {
  map(results, "ks_test") %>%
    bind_rows(.id = "sim_id") %>%
    mutate(
      D = as.numeric(D),
      parameter = factor(parameter, levels = unique(.data$parameter))
    ) %>%
    select(-method2, -l) %>%
    pivot_wider(
      names_from = method1,
      values_from = D
    ) %>%
    ggplot(aes(x = aghq, y = TMB)) +
      geom_point() +
      facet_wrap(~parameter) +
      xlim(0, 0.5) +
      ylim(0, 0.5) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      labs(x = "KS(aghq, tmbstan)", y = "KS(TMB, tmbstan)") +
      theme_minimal()
}

draw_ksplots_l <- function(results) {
  map(results, "ks_test") %>%
    bind_rows(.id = "sim_id") %>%
    mutate(
      l = as.numeric(l),
      parameter = factor(parameter, levels = unique(.data$parameter))
    ) %>%
    select(-D, -method2) %>%
    ggplot(aes(y = l, fill = method1)) +
      geom_boxplot(alpha = 0.7) +
      facet_wrap(~parameter, scales = "free") +
      scale_fill_manual(values = cbpalette) +
      labs(fill = "Method") +
      theme_minimal()
}

draw_rhatplot <- function(results) {
  map(results, "mcmc_monitor") %>%
    map(as.data.frame) %>%
    map(tibble::rownames_to_column, var = "parameter") %>%
    bind_rows(.id = "sim_id") %>%
    ggplot(aes(x = parameter, y = Rhat)) +
    geom_point() +
    geom_hline(yintercept = 1.1, linetype = "dashed") +
    labs(x = "Parameter")  +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90))
}
