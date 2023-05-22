# elgm-inf

[Naomi](https://github.com/mrc-ide/naomi) ([Eaton et al, 2021](https://onlinelibrary.wiley.com/doi/10.1002/jia2.25788)) is a spatial evidence synthesis model used to produce district-level HIV epidemic indicators in sub-Saharan Africa.
Multiple outcomes of interest, including HIV prevalence, HIV incidence and treatment coverage are jointly modelled using both household survey data and routinely reported health system data.
The model is provided as a [tool](https://naomi.unaids.org/) for countries to input their data to and generate estimates using an empirical Bayes Gaussian approximation via the [`TMB`](https://kaskr.github.io/adcomp/_book/Introduction.html) R package.
We propose a new inference method extending adaptive Gauss-Hermite quadrature to deal with >20 hyperparameters, thereby enabling fast and accurate inference for Naomi and other [extended latent Gaussian](https://www.tandfonline.com/doi/full/10.1080/10618600.2022.2099403) models.
Using data from Malawi, our method provides more accurate inferences than `TMB`, and is comparable to Hamiltonian Monte Carlo with the No-U-Turn sampler, but faster to run.
By extending the [`aghq`](https://github.com/awstringer1/aghq) package ([Stringer, 2021](https://arxiv.org/abs/2101.04468)) we facilitate easy, flexible use of our method when provided a [`TMB`](https://kaskr.github.io/adcomp/_book/Introduction.html) C++ template for the model's log-posterior.

![Example district-level Naomi model outputs for adults aged 15-49.](naomi_results.png)

## File structure

The directories of this repository are:

| Directory   | Contains |
|-------------|--------------|
| `docs`      | Hosted files at [https://athowes.github.io/elgm-inf/](https://athowes.github.io/elgm-inf/) |
| `make`      | Scripts used to run the reports. `_make.R` runs everything in order. |
| `src`       | All `orderly` reports. |
| `utils`     | Helper scripts for common development tasks. |

## `orderly`

We use the [`orderly`](https://github.com/vimc/orderly) package ([RESIDE, 2020](https://reside-ic.github.io/)) to simplify the process of doing reproducible research.
After installing [`orderly`](https://github.com/vimc/orderly) (from either CRAN or Github) a report, `example`, may be run by:

```r
orderly::orderly_run(name = "src/example")
```

The results of this run will appear in the `draft/` folder (ignored on Github).
To commit the draft (with associated `id`) to the `archive/` folder (also ignored on Github, and treated as "read only") use:

```r
orderly::orderly_commit(id)
```

Any outputs of this report will then be available to use as dependencies within other reports.
Reports can be pushed to the HIV inference group sharepoint (the remote) using:

```r
orderly::orderly_push_archive("example")
```

Or can be pulled (alongside any dependencies) from the remote using:

```r
orderly_pull_archive("example")
```

Alternatively, just the dependencies can be pulled using `orderly::orderly_pull_dependencies("example")`.

## R package dependencies

This repository is supported by the [`inf.utils`](https://github.com/athowes/inf.utils) package, which can be installed from Github via:

```r
devtools::install_github("athowes/inf.utils")
```

We also use the [`aghq`](https://github.com/awstringer1/aghq) package for inference, which is available from CRAN, though the latest development version can be installed from Github via:

```r
devtools::install_github("athowes/aghq", ref = "adam-dev")
```

## Session information

The `sessionInfo()` used to run this analysis is:

```
R version 4.2.0 (2022-04-22)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS 13.3.1

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] patchwork_1.1.2 stringr_1.5.0   purrr_1.0.1     readr_2.1.3     tidyr_1.2.1     tibble_3.2.1    tidyverse_1.3.1 forcats_0.5.2   ggplot2_3.4.0  
[10] dplyr_1.0.10    rmarkdown_2.18 
```
