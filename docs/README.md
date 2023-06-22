## Notebooks

### Learning and capacity building

* [naomi](https://athowes.github.io/naomi-aghq/naomi.html): following `naomi` vignettes
* [maths](https://athowes.github.io/naomi-aghq/maths.html): mathematical description of the Naomi model
* [aghq](https://athowes.github.io/naomi-aghq/aghq.html): following `aghq` vignettes to explore the package
* [walkthrough](https://athowes.github.io/naomi-aghq/walkthrough.html): stepping through `aghq` code line-by-line
* [`TMB`](https://athowes.github.io/naomi-aghq/tmb.pdf): note about what `TMB` is doing
* [inla-replication](https://athowes.github.io/naomi-aghq/inla-replication.html): replicating the "INLA from scratch" section of [Spatial and Spatio-temporal Bayesian Models with `R-INLA`](https://onlinelibrary.wiley.com/doi/book/10.1002/9781118950203) using R and `TMB`
* [inla-grid](https://athowes.github.io/naomi-aghq/inla-grid.html): illustration of how the hyperparameter posterior marginal is explored in the INLA method
* [epil](https://athowes.github.io/naomi-aghq/epil.html): comparison of Stan, `R-INLA`, `TMB`, `glmmTMB`, `tmbstan` and `aghq` for the epilepsy example from [Rue, Martino and Chopin (2009)](https://rss.onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2008.00700.x)
* [prev-anc-art](https://athowes.github.io/naomi-aghq/prev-anc-art.html): comparison of `tmbstan`, `TMB` and `aghq` for a collection of models based on Joint small-area estimation of HIV prevalence, ART coverage, and HIV incidence (Eaton et al. 2019)

### Grid scale-up

* [scale-grid](https://athowes.github.io/naomi-aghq/scale-grid.html): how can we scale up the number of points in the hyperparameter grid in an intelligent way?
* [astro](https://athowes.github.io/naomi-aghq/astro.html): application of scalable grid to astronomy example from [Bilodeau, Stringer and Tang (2022)](https://www.tandfonline.com/doi/full/10.1080/01621459.2022.2141635)
* [increase-s-k](https://athowes.github.io/naomi-aghq/increase-s-k.html): how does increasing $s$ or $k$ change estimation of the log normalising constant for Naomi?

### Laplace marginals

* [sinla](https://athowes.github.io/naomi-aghq/sinla.html): implementing approximations for the posterior marginal of the latent field, building to the approach of [Wood (2021)](https://academic.oup.com/biomet/article/107/1/223/5572662)

### Posterior comparison approaches

* [posterior-comparison](https://athowes.github.io/naomi-aghq/posterior-comparison.html): exploring methods (Kolmogorov-Smirnov, simulation-based calibration, Pareto smoothed importance sampling, maximum mean discrepancy) for comparison of posterior distributions from approximate Bayesian inference methods

### Results

* [ks](https://athowes.github.io/naomi-aghq/ks.html): comparison of inference methods for the simplified Naomi model using histograms and KS test results
* [exceedance](https://athowes.github.io/naomi-aghq/exceedance.html): case-study of exceedance probabilities: probability to meet the second 90 target, high incidence strata, and amount of unmet treatment need
* [psis](https://athowes.github.io/naomi-aghq/psis.html): comparison of inference methods for the simplified Naomi model using Pareto-smoothed importance sampling
* [mmd](https://athowes.github.io/naomi-aghq/mmd.html): comparison of inference methods for the simplified Naomi model using maximum mean discrepancy
* [mcmc-convergence](https://athowes.github.io/naomi-aghq/mcmc-convergence.html): assessing MCMC (NUTS using `tmbstan`) convergence for the simplified Naomi model
* [model-checks](https://athowes.github.io/naomi-aghq/model-checks.html): checking the fit of the simplified Naomi model to data

<!--

## Experiments

| `TMB` template      | Sample size parameter | Results  |
|:--------------------|:----- |:-----------|
| `model1.cpp`        | 1     | [Plots](https://athowes.github.io/naomi-aghq/model1-plots-m1.pdf) |
| `model1.cpp`        | 10    | [Plots](https://athowes.github.io/naomi-aghq/model1-plots-m10.pdf) |
| `model1.cpp`        | 100   | [Plots](https://athowes.github.io/naomi-aghq/model1-plots-m100.pdf) |
| `model1.cpp`        | 250   | [Plots](https://athowes.github.io/naomi-aghq/model1-plots-m250.pdf) |
| `model1_icar.cpp`   | 1     | [Plots](https://athowes.github.io/naomi-aghq/model1-icar-plots-m1.pdf) |
| `model1_icar.cpp`   | 10    | [Plots](https://athowes.github.io/naomi-aghq/model1-icar-plots-m10.pdf) |
| `model1_icar.cpp`   | 100   | [Plots](https://athowes.github.io/naomi-aghq/model1-icar-plots-m100.pdf) |
| `model1_icar.cpp`   | 250   | [Plots](https://athowes.github.io/naomi-aghq/model1-icar-plots-m250.pdf) |

-->

## Presentations

* [01-04-20-mini](https://athowes.github.io/naomi-aghq/01-04-20-mini.pdf): three month CDT mini-project
* [01-07-21-stats-epi-group](https://athowes.github.io/naomi-aghq/01-07-21-stats-epi-group.pdf): [Machine Learning & Global Health Network](https://mlgh.net/) group meeting
* [15-11-22-seminar](https://athowes.github.io/naomi-aghq/15-11-22-seminar.pdf): Waterloo [SAS Student Seminar Series](https://uwaterloo.ca/statistics-and-actuarial-science/student-seminar-series), and [abstract](https://athowes.github.io/naomi-aghq/seminar.html)
* [18-04-23-lab-group](https://athowes.github.io/naomi-aghq/18-04-23-lab-group.pdf): HIV inference group lab meeting
* [19-05-23-turing](https://athowes.github.io/naomi-aghq/19-05-23-turing.pdf): [Alan Turing Insitute PhD Student Presentation and Networking](https://www.turing.ac.uk/events/phd-student-presentation-and-networking-may-2023)
* [21-06-23-mlgh](https://athowes.github.io/naomi-aghq/21-06-23-mlgh.pdf): [Machine Learning & Global Health Network](https://mlgh.net/) group meeting

## Posters

* [bayescomp-poster](https://athowes.github.io/naomi-aghq/bayescomp-poster.pdf): for BayesComp 2023
* [bioinference-poster](https://athowes.github.io/naomi-aghq/bioinference-poster.pdf): for BioInference 2023

## Manuscript

* [paper](https://athowes.github.io/naomi-aghq/paper.pdf): work-in-progress write-up
* [appendix](https://athowes.github.io/naomi-aghq/appendix.pdf): additional material for write-up

## Misc

* [bayescomp](https://athowes.github.io/naomi-aghq/bayescomp.html): abstract for BayesComp 2023
* [bioinference](https://athowes.github.io/naomi-aghq/bioinference.html): abstract for BioInference 2023
