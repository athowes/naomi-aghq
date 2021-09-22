draw_boxplots <- function(results, methods) {
  map(results, "comparison_results") %>%
    bind_rows(.id = "sim_id") %>%
    pivot_longer(cols = all_of(methods), names_to = "method", values_to = "value") %>%
    ggplot(aes(x = method, y = value, fill = method)) +
      geom_boxplot() +
      facet_wrap(parameter ~ type, scales = "free") +
      scale_fill_manual(values = cbpalette) +
      labs(x = "", y = "Posterior mean", fill = "Inference method")
}
