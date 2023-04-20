#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("check_tmb-output")
# setwd("src/check_tmb-output")

tmb <- readRDS("depends/tmb.rds")

qq <- function(data) {
  qqnorm(data)
  qqline(data)
}

pdf("qq-checks.pdf", h = 4, w = 6.25)

qq(tmb$fit$sample$beta_alpha[1, ])
qq(tmb$fit$sample$beta_alpha[2, ])

qq(tmb$fit$sample$u_rho_xs[1, ])
qq(tmb$fit$sample$u_rho_xs[2, ])

qq(tmb$fit$sample$rho_t1_out[1, ])
qq(tmb$fit$sample$rho_t1_out[2, ])

qq(tmb$fit$sample$alpha_t1_out[1, ])
qq(tmb$fit$sample$alpha_t1_out[2, ])

dev.off()
