#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("random-mixture")
# setwd("src/random-mixture")

compile("random-mixture.cpp", flags = "-w")
dyn.load(dynlib("random-mixture"))
