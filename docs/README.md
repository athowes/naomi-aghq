## Notebooks

* [aghq](https://athowes.github.io/elgm-inf/aghq.html): following `aghq` vignettes to explore the package
* [walkthrough](https://athowes.github.io/elgm-inf/walkthrough.html): stepping through `aghq` code line-by-line
* [epil](https://athowes.github.io/elgm-inf/epil.html): comparison of Stan, `R-INLA`, `TMB`, `glmmTMB`, `tmbstan` and `aghq` for the epilepsy example from [Rue, Martino and Chopin (2009)](https://rss.onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2008.00700.x)
* [inla-grid](https://athowes.github.io/elgm-inf/inla-grid.html): illustration of how the hyperparameter posterior marginal is explored in the INLA method
* [inla-replication](https://athowes.github.io/elgm-inf/inla-replication.html): replicating the "INLA from scratch" section of [Spatial and Spatio-temporal Bayesian Models with `R-INLA`](https://onlinelibrary.wiley.com/doi/book/10.1002/9781118950203) using R and `TMB`
* [sinla](https://athowes.github.io/elgm-inf/sinla.html): implementing approximations for the posterior marginal of the latent field, building to the approach of [Wood (2021)](https://academic.oup.com/biomet/article/107/1/223/5572662)
* [prev-anc-art](https://athowes.github.io/elgm-inf/prev-anc-art.html): comparison of `tmbstan`, `TMB` and `aghq` for a collection of models based on Joint small-area estimation of HIV prevalence, ART coverage, and HIV incidence (Eaton et al. 2019)
* [naomi](https://athowes.github.io/elgm-inf/naomi.html): following `naomi` vignettes
* [maths](https://athowes.github.io/elgm-inf/maths.html): mathematical description of the Naomi model
* [posterior-comparison](https://athowes.github.io/elgm-inf/posterior-comparison.html): exploring methods (Kolmogorov-Smirnov, simulation-based calibration, Pareto smoothed importance sampling, maximum mean discrepancy) for comparison of posterior distributions from approximate Bayesian inference methods
* [mcmc-convergence](https://athowes.github.io/elgm-inf/mcmc-convergence.html): assessing MCMC (NUTS using `tmbstan`) convergence for the simplified Naomi model
* [comparison](https://athowes.github.io/elgm-inf/comparison.html): comparison of inference methods for the simplified Naomi model
* [model-checks](https://athowes.github.io/elgm-inf/model-checks.html): checking the fit of the simplified Naomi model to data
* [aghq-grid-scale-up](https://athowes.github.io/elgm-inf/aghq-grid-scale-up.html): how can we scale up the number of points in the hyperparameter grid in an intelligent way?

## Experiments

| `TMB` template      | Sample size parameter | Results  |
|:--------------------|:----- |:-----------|
| `model1.cpp`        | 1     | [Plots](https://athowes.github.io/elgm-inf/model1-plots-m1.pdf) |
| `model1.cpp`        | 10    | [Plots](https://athowes.github.io/elgm-inf/model1-plots-m10.pdf) |
| `model1.cpp`        | 100   | [Plots](https://athowes.github.io/elgm-inf/model1-plots-m100.pdf) |
| `model1.cpp`        | 250   | [Plots](https://athowes.github.io/elgm-inf/model1-plots-m250.pdf) |
| `model1_icar.cpp`   | 1     | [Plots](https://athowes.github.io/elgm-inf/model1-icar-plots-m1.pdf) |
| `model1_icar.cpp`   | 10    | [Plots](https://athowes.github.io/elgm-inf/model1-icar-plots-m10.pdf) |
| `model1_icar.cpp`   | 100   | [Plots](https://athowes.github.io/elgm-inf/model1-icar-plots-m100.pdf) |
| `model1_icar.cpp`   | 250   | [Plots](https://athowes.github.io/elgm-inf/model1-icar-plots-m250.pdf) |

## Presentations

* [01-04-20-mini](https://athowes.github.io/elgm-inf/01-04-20-mini.pdf): three month CDT mini-project
* [01-07-21-stats-epi-group](https://athowes.github.io/elgm-inf/01-07-21-stats-epi-group.pdf): [Machine Learning & Global Health Network](https://mlgh.net/) group meeting
* [15-11-22-seminar](https://athowes.github.io/elgm-inf/15-11-22-seminar.pdf): Waterloo [SAS Student Seminar Series](https://uwaterloo.ca/statistics-and-actuarial-science/student-seminar-series)
  * [seminar](https://athowes.github.io/elgm-inf/seminar.html): abstract for Waterloo SAS Student Seminar Series
* [xx-12-22-explainer](https://athowes.github.io/elgm-inf/xx-12-22-explainer.pdf): HIV inference group lab meeting

## Posters

* [bayescomp-poster](https://athowes.github.io/elgm-inf/bayescomp-poster.pdf): for BayesComp 2023

## Manuscript

* [paper](https://athowes.github.io/elgm-inf/paper.pdf): work-in-progress write-up
* [appendix](https://athowes.github.io/elgm-inf/appendix.pdf): additional material for write-up

## Misc

* [bayescomp](https://athowes.github.io/elgm-inf/bayescomp.html): abstract for BayesComp 2023
* [`TMB`](https://athowes.github.io/elgm-inf/tmb.pdf): note about what `TMB` is doing
