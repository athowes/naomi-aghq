# elgm-inf

Given a user template `model.cpp` we compare approaches to Bayesian inference including empirical Bayes via [`TMB`](https://kaskr.github.io/adcomp/Introduction.html), adaptive Gaussian-Hermite quadrature via [`aghq`](https://github.com/awstringer1/aghq) and Hamiltonian Monte Carlo via [`tmbstan`](https://github.com/kaskr/tmbstan).
The motivation for this work is to improve the inference in the [Naomi model for subnational HIV estimates](https://github.com/mrc-ide/naomi).
This may be done by developing an alternative implementation of the integrated nested Laplace approximation (to [`R-INLA`](https://www.r-inla.org/)) by extending the `aghq` package.

## R package dependencies

This repository is supported by the [`inf.utils`](https://github.com/athowes/inf.utils) package, which can be installed from Github via:

```r
devtools::install_github("athowes/inf.utils")
```
