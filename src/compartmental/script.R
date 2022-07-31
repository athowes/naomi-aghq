#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("compartmental")
# setwd("src/compartmental")

compile("compartmental.cpp", flags = "-w")
dyn.load(dynlib("compartmental"))

b <- 0.4
g <- 0.04
I_0 <- 1e-4
max_T <- 100
dt <- 1/10

in_dat <- list(
  N_t = max_T/dt,
  dt = dt
)

in_par <- list(
  b = b,
  g = g,
  I_0 = I_0
)

obj <- MakeADFun(in_dat, in_par, DLL = "compartmental")
matplot(obj$report()$Z_out, type = "l")
