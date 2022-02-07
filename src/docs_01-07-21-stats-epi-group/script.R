#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("docs_01-07-21-stats-epi-group")
# setwd("src/docs_01-07-21-stats-epi-group")

sf_lightgrey <- "#E6E6E6"
lightgrey <- "#D3D3D3"

lightblue <- "#56B4E9"
lightgreen <- "#009E73"
lightpink <- "#CC79A7"
light_palette <- c(lightblue, lightgreen, lightpink)

midblue <- "#3D9BD0"
midgreen <- "#00855A"
midpink <- "#B3608E"
mid_palette <- c(midblue, midgreen, midpink)

darkblue <- "#004E83"
darkgreen <- "#00380D"
darkpink <- "#802D5B"
dark_palette <- c(darkblue, darkgreen, darkpink)

cbpalette <- multi.utils::cbpalette()

rmarkdown::render("01-07-21-stats-epi-group.Rmd")
