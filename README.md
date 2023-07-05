# naomi-aghq

Code for the manuscript Howes, Stringer, Flaxman and Eaton "Fast approximate Bayesian inference of HIV indicators using PCA adaptive Gauss-Hermite quadrature" (in preparation).

[Naomi](https://github.com/mrc-ide/naomi) ([Eaton et al, 2021](https://onlinelibrary.wiley.com/doi/10.1002/jia2.25788)) is a spatial evidence synthesis model used to produce district-level HIV epidemic indicators in sub-Saharan Africa.
Multiple outcomes of interest, including HIV prevalence, HIV incidence and treatment coverage are jointly modelled using both household survey data and routinely reported health system data.
The model is provided as a [tool](https://naomi.unaids.org/) for countries to input their data to and generate estimates during a yearly process supported by UNAIDS.
Currently, inference is conducted using empirical Bayes and a Gaussian approximation via the [`TMB`](https://kaskr.github.io/adcomp/_book/Introduction.html) R package.
We propose a new inference method extending adaptive Gauss-Hermite quadrature to deal with >20 hyperparameters, enabling fast and accurate inference for Naomi and other [extended latent Gaussian](https://www.tandfonline.com/doi/full/10.1080/10618600.2022.2099403) models.
Using data from Malawi, our method improves the accuracy of inferences across a range of metrics, while being substantially faster to run than Hamiltonian Monte Carlo with the No-U-Turn sampler.
Our implementation is based on the existing [`TMB`](https://kaskr.github.io/adcomp/_book/Introduction.html) C++ template for the model's log-posterior, and is compatible with any model with such a template.

![Example district-level Naomi model outputs for adults aged 15-49.](naomi_results.png)

## C++ template for the log-posterior

The `TMB` template for the simplified Naomi model is available [here](https://github.com/athowes/naomi-aghq/blob/master/src/naomi-simple_fit/naomi_simple.cpp).

## File structure

The directories of this repository are:

| Directory | Contains |
|-----------|----------|
| `docs`    | Hosted files at [https://athowes.github.io/naomi-aghq/](https://athowes.github.io/naomi-aghq/) |
| `make`    | Scripts used to run the reports. `_make.R` runs everything in order |
| `misc`    | Ideas for further work and other documents that don't need to be in `src` |
| `src`     | All `orderly` reports |
| `utils`   | Helper scripts for common development tasks |

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

## Reports

The reports within `src` are as follows:

| Report | Description |
|--------|-------------|
| `astro` | Analysis of an astronomy example |
| `check_hyper-marginals` | Visualise hyperparameter posteriors fitted with NUTS for simplified Naomi, as well as positions of PCA-AGHQ nodes |
| `check_pca-aghq` | Unit tests for the PCA-AGHQ procedure working as expected |
| `check_sd-estimation` | Analysing odd behaviour for simplified Naomi whereby joint Gaussian approximation better estimates the standard deviation than conditionals |
| `check_tmb-aghq-k1` | Assessment that TMB and `aghq` run with `k = 1` produce essentially the same inferences for simplified Naomi, in contrast to inferences from `aghq` with any grid |
| `check_tmb-output` | Checking the the TMB output for simplified Naomi looks Gaussian, as it should be |
| `dev_hyper-sampling` | Method development for sampling from the hyperparameters within TMB or AGHQ |
| `dev_scale-grid` | Method development for hyperparameter grids over many dimensional spaces |
| `dev_sinla` | Method development for Laplace and reduced cost Laplace marginals |
| `docs_01-04-20-mini` | Presentation for StatML CDT |
| `docs_01-07-21-stats-epi-group` | Presentation for Flaxman, Bhatt groups |
| `docs_15-11-22-seminar` | Presentation for Waterloo SAS |
| `docs_18-04-23-lab-group` | Presentation for HIV Inference Group |
| `docs_19-05-23-turing` | Presentation for Alan Turing Institute |
| `docs_21-06-23-mlgh` | Presentation for Machine Learning and Global Health Network |
| `docs_bayescomp-poster`| Poster for BayesComp 2023 |
| `docs_bioinference-poster`| Poster for BioInference 2023 |
| `docs_model-structure` | Analysis of inputs to simplified Naomi outputs |
| `docs_paper` | Paper and appendix |
| `epil` | Analysis of epilepsy example |
| `example_inla-grid` | Demonstration of grid construction method used in `R-INLA` |
| `example_inla-replication` | Step-by-step version of INLA method in R and `TMB` |
| `example_naomi` | Follow the `naomi` package vignette |
| `explore_aghq` | Walkthrough various `aghq` package functions |
| `explore_posterior-comparison` | Explore methods for comparing the accuracy of computed posterior distributions |
| `naomi-simple_contraction` | Comparison of prior and posterior standard deviations across inference methods |
| `naomi-simple_exceedance` | Case-study on computation of exceedance probabilities for second 90 target and high incidence |
| `naomi-simple_fit` | Fit the simplified Naomi model using a range of inference methods |
| `naomi-simple_increase-s-k` | Is the log normalising constant better estimated as $s$ and $k$ are increased in PCA-AGHQ? |
| `naomi-simple_ks` | Comparison of inference methods using Kolmogorov-Smirnov tests on marginals |
| `naomi-simple_mcmc` | Diagnostics for MCMC convergence |
| `naomi-simple_mmd` | Comparison of inference methods using maximum mean discrepancy on joint distributions |
| `naomi-simple_model-checks` | Checking the model fit to data |
| `naomi-simple_point-estimates` | Comparison of inference methods for estimating the mean and standard deviation |
| `naomi-simple_psis` | Comparison of inference methods using Pareto-smoothed importance sampling on joint distributions |
| `plot-tikz_algorithm-flowchart` | TikZ diagram of proposed algorithm |
| `plot-tikz_simplified-naomi` | TikZ diagram of the simplified Naomi directed acyclic graph |
| `prev-anc-art_model0` | Fit prevalence, ANC, ART model version 0 |
| `prev-anc-art_model1` | Fit prevalence, ANC, ART model version 1 |
| `prev-anc-art_model2` | Fit prevalence, ANC, ART model version 2 |
| `prev-anc-art_model3` | Fit prevalence, ANC, ART model version 3 |
| `prev-anc-art_model4` | Fit prevalence, ANC, ART model version 4 |
| `prev-anc-art_results` | Analyse results of prevalence, ANC, ART model |
| `prev-anc-art_sim` | Simulate prevalence, ANC, ART data |
