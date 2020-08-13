
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rctSurv

<!-- badges: start -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
<!-- badges: end -->

> Efficient Estimators for Survival Analysis in RCTs Without
> Proportional Hazards

# Installation

The development version can be installed from
[GitHub](https://github.com) with:

``` r
devtools::install_github("nt-williams/rctSurv")
```

# Scope

`rctSurv` implements efficient estimators for the restricted mean
survival time (RMST) and survival probability in randomized controlled
trials (RCT) without the proportional hazards assumption. The provided
estimators are non-parametric and can incorporate flexible,
data-adaptive estimation (provided using the SuperLearner from the `sl3`
package) for nuisance parameters while remaining root-n consistent.

# Example

Using the `veteran` dataset provided by the `survival` package…

``` r
library(rctSurv)

veteran <- survival::veteran
veteran$trt <- veteran$trt - 1
veteran$celltype <- factor(veteran$celltype)
veteran$id <- 1:nrow(veteran)
```

Defining data structure…

``` r
trt <- "trt"
time <- "time"
cens <- "status"
covar <- c("celltype", "karno", "diagtime", "age", "prior")
id <- "id"
```

Setting up the SuperLearner for estimation…

``` r
lrnrs <- sl3::make_learner_stack(sl3::Lrnr_glm_fast, 
                                 sl3::Lrnr_ranger)
```

Estimate nuisance parameters and establish metadata…

``` r
meta <- metadata(veteran, trt, cens, covar, time, id,
                 coarsen = 30, estimator = "tmle", 
                 lrnrs_trt = lrnrs, lrnrs_cens = lrnrs, lrnrs_hzrd = lrnrs)
meta
#> rctSurv metadata
#> 
#> ● Estimate RMST with `rmst()`
#> ● Estimate survival probability with `survprob()`
#> ● Inspect SuperLearner weights with `get_weights()`
#> 
#>          Estimator: tmle
#>   Treatment status: trt
#>   Censoring status: status
#>     Adjustment set: celltype, karno, diagtime, age, and prior
#> Max coarsened time: 34
```

We can compute restricted mean survival time…

``` r
rmst(meta)
#> RMST Estimator: tmle
#> 
#>            Confidence level: 95%
#>      Multiplier Bootstrap C: 2.401891 
#>  Test of no effect, p-value:
#>       First 6 time horizons:
#> 
#>   horizon Treatment Control Theta Point-wise 95% CI  Uniform 95% CI
#> 1       2      1.69    1.73 -0.03   (-0.17 to 0.10) (-0.20 to 0.13)
#> 2       3      2.19    2.29 -0.11   (-0.35 to 0.14) (-0.41 to 0.19)
#> 3       4      2.59    2.82 -0.23   (-0.59 to 0.13) (-0.67 to 0.21)
#> 4       5      2.87    3.22 -0.34   (-0.80 to 0.12) (-0.91 to 0.23)
#> 5       6      3.12    3.52 -0.40   (-0.97 to 0.17) (-1.10 to 0.30)
#> 6       7      3.35    3.74 -0.39   (-1.05 to 0.27) (-1.20 to 0.42)
```

We can also compute survival probabilities…

``` r
survprob(meta)
#> Survival Probability Estimator: tmle
#> 
#>            Confidence level: 95%
#>      Multiplier Bootstrap C: 2.699171 
#>  Test of no effect, p-value:
#>       First 6 time horizons:
#> 
#>   horizon Treatment Control Theta Point-wise 95% CI  Uniform 95% CI
#> 1       2      0.69    0.73 -0.03   (-0.17 to 0.10) (-0.22 to 0.15)
#> 2       3      0.49    0.57 -0.08   (-0.22 to 0.07) (-0.28 to 0.13)
#> 3       4      0.40    0.52 -0.12   (-0.26 to 0.03) (-0.32 to 0.08)
#> 4       5      0.29    0.40 -0.11   (-0.27 to 0.04) (-0.32 to 0.10)
#> 5       6      0.25    0.31 -0.06   (-0.21 to 0.09) (-0.27 to 0.15)
#> 6       7      0.23    0.22  0.01   (-0.12 to 0.14) (-0.17 to 0.19)
#> Access all estimates with `all_estimates()`
```

# Citation

# References

Díaz, I., E. Colantuoni, D. F. Hanley, and M. Rosenblum (2019). Improved
precision in the analysis of randomized trials with survival outcomes,
without assuming proportional hazards. Lifetime Data Analysis 25 (3),
439–468.
