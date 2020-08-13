
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
```

Setting up the SuperLearner for estimation…

``` r
lrnrs <- sl3::make_learner_stack(sl3::Lrnr_glm_fast, 
                                 sl3::Lrnr_ranger)
```

Estimate nuisance parameters and establish metadata…

``` r
surv <- survrct(Surv(time, status) ~ trt + celltype + karno + diagtime + age + prior, 
                target = "trt", data = veteran, coarsen = 7, estimator = "tmle", 
                lrnrs_trt = lrnrs, lrnrs_cens = lrnrs, lrnrs_hzrd = lrnrs)
#> survrct metadata
#> 
#> Surv(time, status) ~ trt + celltype + karno + diagtime + age + 
#>     prior
#> 
#> ● Estimate RMST with `rmst()`
#> ● Estimate survival probability with `survprob()`
#> ● Inspect SuperLearner weights with `get_weights()`
#> 
#>          Estimator: tmle
#>    Target variable: trt
#>   Status Indicator: status
#>     Adjustment set: celltype, karno, diagtime, age, and prior
#> Max coarsened time: 143
```

We can compute restricted mean survival time…

``` r
rmst(surv)
#> RMST Estimator: tmle
#> 
#>            Confidence level: 95%
#>      Multiplier Bootstrap C: 2.530097 
#>  Test of no effect, p-value:
#>       First 6 time horizons:
#> 
#>   horizon Treatment Control Theta Point-wise 95% CI  Uniform 95% CI
#> 1       2      1.96    1.97 -0.01   (-0.08 to 0.05) (-0.09 to 0.07)
#> 2       3      2.84    2.80  0.04   (-0.11 to 0.20) (-0.16 to 0.24)
#> 3       4      3.65    3.58  0.07   (-0.19 to 0.33) (-0.26 to 0.40)
#> 4       5      4.36    4.30  0.06   (-0.30 to 0.42) (-0.41 to 0.53)
#> 5       6      5.02    5.00  0.02   (-0.46 to 0.49) (-0.60 to 0.63)
#> 6       7      5.66    5.68 -0.02   (-0.62 to 0.57) (-0.79 to 0.74)
#> Access all estimates with `all_estimates()`
```

We can also compute survival probabilities…

``` r
survprob(surv)
#> Survival Probability Estimator: tmle
#> 
#>            Confidence level: 95%
#>      Multiplier Bootstrap C: 2.877887 
#>  Test of no effect, p-value:
#>       First 6 time horizons:
#> 
#>   horizon Treatment Control Theta Point-wise 95% CI  Uniform 95% CI
#> 1       2      0.96    0.97 -0.01   (-0.08 to 0.05) (-0.11 to 0.08)
#> 2       3      0.89    0.83  0.05   (-0.07 to 0.17) (-0.12 to 0.23)
#> 3       4      0.80    0.77  0.03   (-0.10 to 0.15) (-0.16 to 0.21)
#> 4       5      0.71    0.73 -0.01   (-0.15 to 0.12) (-0.21 to 0.18)
#> 5       6      0.66    0.70 -0.04   (-0.18 to 0.10) (-0.25 to 0.16)
#> 6       7      0.64    0.68 -0.04   (-0.18 to 0.10) (-0.24 to 0.16)
#> Access all estimates with `all_estimates()`
```

# Citation

# References

Díaz, I., E. Colantuoni, D. F. Hanley, and M. Rosenblum (2019). Improved
precision in the analysis of randomized trials with survival outcomes,
without assuming proportional hazards. Lifetime Data Analysis 25 (3),
439–468.