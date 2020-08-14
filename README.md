
<!-- README.md is generated from README.Rmd. Please edit that file -->

# survrct

<!-- badges: start -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Travis build
status](https://travis-ci.com/nt-williams/survrct.svg?branch=master)](https://travis-ci.com/nt-williams/survrct)
[![Codecov test
coverage](https://codecov.io/gh/nt-williams/rctSurv/branch/master/graph/badge.svg)](https://codecov.io/gh/nt-williams/survrct?branch=master)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

> Efficient Estimators for Survival Analysis in RCTs Without
> Proportional Hazards

Nick Williams and Ivan Diaz

-----

# Installation

The development version can be installed from
[GitHub](https://github.com) with:

``` r
devtools::install_github("nt-williams/survrct")
```

# Scope

`survrct` implements efficient estimators for the restricted mean
survival time (RMST) and survival probability in randomized controlled
trials (RCT) without the proportional hazards assumption as introduced
in Diaz et al. (2019). Prognostic baseline variables should be
incorporated to obtain equal or better asymptotic precision compared to
unadjusted estimators.

Under random censoring, the primary estimator (TMLE) is doubly robust–it
is consistent if either the outcome or censoring model is correctly
specified. This is an advantage over Cox proportional hazards which
relies on assumptions of only the outcome model. The primary estimator
is non-parametric and can incorporate flexible, data-adaptive estimation
(provided using the SuperLearner from the
[`sl3`](https://github.com/tlverse/sl3) package) for nuisance parameters
while remaining root-n consistent.

Estimators are implemented using a formula interface based on that of
the [`survival`](https://CRAN.R-project.org/package=survival) package
for familiar users.

# Example

Using the `veteran` dataset provided by the
[`survival`](https://CRAN.R-project.org/package=survival) package:

``` r
library(survrct)
library(future)

veteran <- survival::veteran
veteran$trt <- veteran$trt - 1
veteran$celltype <- factor(veteran$celltype)
head(veteran)
#>   trt celltype time status karno diagtime age prior
#> 1   0 squamous   72      1    60        7  69     0
#> 2   0 squamous  411      1    70        5  64    10
#> 3   0 squamous  228      1    60        3  38     0
#> 4   0 squamous  126      1    60        9  63    10
#> 5   0 squamous  118      1    70       11  65    10
#> 6   0 squamous   10      1    20        5  49     0
```

We can use the SuperLearner provided by the
[`sl3`](https://github.com/tlverse/sl3) package for estimation of the
outcome and censoring models. Alternatively, if no SuperLearner
specification is provided, estimation will be performed using saturated
GLMs.

To allow for estimation of multiple estimands without having to
re-estimate nuisance parameters we first create a Survival metadata
object using the `survrct()` function. We specify the model
paramaterization using a typical `R` formula with `Surv()` (based on the
[`survival`](https://CRAN.R-project.org/package=survival) package)
specifying the left-hand side of the formula.

We also specify our target variable of interest, an optional time
coarsener, and the estimator we would like to use.

``` r
plan(multisession)
surv <- survrct(Surv(time, status) ~ trt + celltype + karno + diagtime + age + prior, 
                target = "trt", data = veteran, coarsen = 7, estimator = "tmle")
surv
#> survrct metadata
#> 
#> Surv(time, status) ~ trt + celltype + karno + diagtime + age + 
#>     prior
#> 
#> ● Estimate RMST with `rmst()`
#> ● Estimate survival probability with `survprob()`
#> ● Inspect nuisance parameter models with `get_fits()`
#> 
#>          Estimator: tmle
#>             Engine: sl3
#>    Target variable: trt
#>   Status Indicator: status
#>     Adjustment set: celltype, karno, diagtime, age, and prior
#> Max coarsened time: 143
```

Using the metadata from the previous step we can now estimate the
restricted mean survival time. Two confidence bands are returned: 95%
point-wise intervals as well as 95% uniform confidence bands based on
the multiplier-bootstrap from Kennedy (2019).

``` r
rmst(surv)
#> RMST Estimator: tmle
#> 
#>            Confidence level: 95%
#>      Multiplier Bootstrap C: 2.512656 
#>  Test of no effect, p-value:
#>           First 6 estimates:
#> 
#>   horizon Treatment Control Theta Point-wise 95% CI  Uniform 95% CI
#> 1       2      1.96    1.97 -0.01   (-0.08 to 0.05) (-0.09 to 0.07)
#> 2       3      2.84    2.80  0.04   (-0.11 to 0.20) (-0.16 to 0.24)
#> 3       4      3.65    3.58  0.07   (-0.19 to 0.33) (-0.26 to 0.40)
#> 4       5      4.36    4.30  0.06   (-0.31 to 0.42) (-0.41 to 0.52)
#> 5       6      5.02    5.00  0.02   (-0.46 to 0.49) (-0.60 to 0.63)
#> 6       7      5.66    5.68 -0.02   (-0.62 to 0.57) (-0.79 to 0.74)
#> Access all estimates with `all_estimates()`
```

We can also estimate survival probabilities.

``` r
survprob(surv)
#> Survival Probability Estimator: tmle
#> 
#>            Confidence level: 95%
#>      Multiplier Bootstrap C: 2.862583 
#>  Test of no effect, p-value:
#>           First 6 estimates:
#> 
#>   horizon Treatment Control Theta Point-wise 95% CI  Uniform 95% CI
#> 1       1      0.96    0.97 -0.01   (-0.08 to 0.05) (-0.11 to 0.08)
#> 2       2      0.89    0.83  0.05   (-0.07 to 0.17) (-0.12 to 0.23)
#> 3       3      0.80    0.77  0.03   (-0.10 to 0.15) (-0.16 to 0.21)
#> 4       4      0.71    0.73 -0.01   (-0.15 to 0.12) (-0.21 to 0.18)
#> 5       5      0.66    0.70 -0.04   (-0.18 to 0.10) (-0.25 to 0.16)
#> 6       6      0.64    0.68 -0.04   (-0.18 to 0.10) (-0.24 to 0.16)
#> Access all estimates with `all_estimates()`
```

# Citation

Please include the following citations after use:

# References

Díaz, I., E. Colantuoni, D. F. Hanley, and M. Rosenblum (2019). Improved
precision in the analysis of randomized trials with survival outcomes,
without assuming proportional hazards. Lifetime Data Analysis 25 (3),
439–468.

Edward H. Kennedy (2019) Nonparametric Causal Effects Based on
Incremental Propensity Score Interventions, Journal of the American
Statistical Association, 114:526, 645-656, DOI:
10.1080/01621459.2017.1422737
