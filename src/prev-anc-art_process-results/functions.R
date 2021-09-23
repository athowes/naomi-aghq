draw_boxplots <- function(results, methods) {
  map(results, "comparison_results") %>%
    bind_rows(.id = "sim_id") %>%
    pivot_longer(cols = all_of(methods), names_to = "method", values_to = "value") %>%
    left_join(true_values, by = "parameter") %>%
    mutate(true_value = ifelse(type == "SD", NA, true_value)) %>%
    ggplot(aes(x = method, y = value, fill = method)) +
      geom_boxplot() +
      geom_hline(aes(yintercept = true_value), linetype = "dashed") +
      facet_wrap(parameter ~ type, scales = "free", ncol = 2) +
      scale_fill_manual(values = cbpalette) +
      labs(x = "", y = "Estimate", fill = "Inference method")
}

draw_scatterplots <- function(results, methods) {
  df <- map(results, "comparison_results") %>%
    bind_rows(.id = "sim_id") %>%
    pivot_longer(cols = all_of(methods), names_to = "method", values_to = "value") %>%
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
    labs(x = "Gold standard (tmbstan)", y = "Value", col = "Inference method")
}
