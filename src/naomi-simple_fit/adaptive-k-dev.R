#' Given TMB fit chose k for each dimension
hypers <- names(fit$par)

#' Should be 24
length(hypers)

means <- vector(mode = "numeric", length = length(hypers))
sds <- vector(mode = "numeric", length = length(hypers))

for(i in seq_along(hypers)) {
  means[i] <- mean(fit$sample[[hypers[i]]])
  sds[i] <- sd(fit$sample[[hypers[i]]])
}

#' TODO Turn this into an interval plot with parameter names along one side
#' Maybe needs a coord_flip()
data.frame(name = hypers, mean = means, sd = sds) %>%
  ggplot(aes(x = mean, y = sd)) +
  geom_point()
