
<!-- README.md is generated from README.Rmd. Please edit that file -->

# adjrct

<!-- badges: start -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R build
status](https://github.com/nt-williams/rctSurv/workflows/R-CMD-check/badge.svg)](https://github.com/nt-williams/rctSurv/actions)

<!-- badges: end -->

> Efficient Estimators for Survival and Ordinal Outcomes in RCTs Without
> Proportional Hazards and Odds Assumptions with Variable Selection

Nick Williams and Iván Díaz

------------------------------------------------------------------------

# Installation

The development version can be installed from
[GitHub](https://github.com) with:

``` r
devtools::install_github("nt-williams/adjrct")
```

# Scope

**adjrct** implements efficient estimators for the restricted mean
survival time (RMST) and survival probability for time-to-event outcomes
(Díaz et al., 2019), and the average log odds ratio (Díaz et al., 2016)
and Mann-Whitney estimand (Vermeulen et al., 2014) for ordinal outcomes
in randomized controlled trials (RCT) without the proportional hazards
or odds assumptions. Prognostic baseline variables should be
incorporated to obtain equal or better asymptotic precision compared to
un-adjusted estimators. Under random censoring, the primary estimator
(TMLE) is doubly robust–it is consistent if either the outcome or
censoring model is correctly specified. For survival outcomes,
estimators are implemented using a formula interface based on that of
the [**survival**](https://CRAN.R-project.org/package=survival) package
for familiar users.

# Example

### Survival outcome

To allow for estimation of multiple estimands without having to
re-estimate nuisance parameters we first create a Survival metadata
object using the `survrct()` function. We specify the model
parameterization using a typical `R` formula with `Surv()` (based on the
[**survival**](https://CRAN.R-project.org/package=survival) package)
specifying the left-hand side of the formula. We also specify a formula
for the propensity.

``` r
library(adjrct)

data("c19.tte")
surv <- survrct(Surv(days, event) ~ A + age + sex + bmi + dyspnea, 
                A ~ 1, data = c19.tte)
surv
#> survrct metadata
#> 
#> Surv(days, event) ~ A + age + sex + bmi + dyspnea
#> 
#> • Estimate RMST with `rmst()`
#> • Estimate survival probability with `survprob()`
#> • Inspect nuisance parameter models with `get_fits()`
#> 
#>          Estimator: tmle
#>    Target variable: A
#>   Status Indicator: event
#>     Adjustment set: age, sex, bmi, and dyspnea
#> Max coarsened time: 15
```

Using the metadata from the previous step we can now estimate the
restricted mean survival time and survival probability for a single or
multiple time horizons. If multiple times are evaluated, two confidence
bands are returned: 95% point-wise intervals as well as 95% uniform
confidence bands based on the multiplier-bootstrap from Kennedy (2019).

``` r
rmst(surv, 14)
#> RMST Estimator: tmle
#> 
#> Marginal RMST: E(min[T, 14] | A = a)
#> Treatment Arm
#>       Estimate: 12.43
#>     Std. error: 0.12
#>         95% CI: (12.2, 12.66)
#> Control Arm
#>       Estimate: 11.37
#>     Std. error: 0.17
#>         95% CI: (11.04, 11.7)
#> 
#> Treatment Effect: E(min[T, 14] | A = 1) - E(min[T, 14] | A = 0)
#> Additive effect
#>       Estimate: 1.06
#>     Std. error: 0.2
#>         95% CI: (0.66, 1.46)
survprob(surv, 14)
#> Survival Probability Estimator: tmle
#> 
#> Marginal Survival Probability: Pr(T > 14 | A = a)
#> Treatment Arm
#>       Estimate: 0.76
#>     Std. error: 0.02
#>         95% CI: (0.73, 0.79)
#> Control Arm
#>       Estimate: 0.73
#>     Std. error: 0.02
#>         95% CI: (0.7, 0.76)
#> 
#> Treatment Effect: Pr(T > 14 | A = 1) - Pr(T > 14 | A = 0)
#> Additive effect
#>       Estimate: 0.03
#>     Std. error: 0.02
#>         95% CI: (-0.02, 0.07)
```

### Ordinal outcome

We can similarly create an Ordinal metadata object using the
`ordinalrct()` function.

``` r
data("c19.ordinal")

ord <- ordinalrct(state_ordinal ~ A + age + dyspnea + sex, 
                  A ~ 1, data = c19.ordinal)
ord
#> ordinalrct metadata
#> 
#> Outcome regression: state_ordinal ~ A + age + dyspnea + sex
#>         Propensity: A ~ 1
#> 
#> • Estimate log odds ratio with `log_or()`
#> • Estimate Mann-Whitney with `mannwhitney()`
#> • Estimate with `cdf()`
#> • Estimate with `pmf()`
#> • Inspect nuisance parameter models with `get_fits()`
#> 
#>          Estimator: tmle
#>    Target variable: A
#>   Outcome variable: state_ordinal
```

The average log odds ratio, Mann-Whitney statistic, CDF, and PMF can
then be estimated using the metadata.

``` r
log_or(ord)
#> Log OR Estimator: tmle
#> 
#> Arm-specific log odds:
#> Treatment Arm
#>       Estimate: 1.62
#>     Std. error: 0.09
#>         95% CI: (1.44, 1.8)
#> Control Arm
#>       Estimate: 0.87
#>     Std. error: 0.07
#>         95% CI: (0.73, 1.02)
#> 
#> Average log odds ratio:
#>       Estimate: 0.75
#>     Std. error: 0.11
#>         95% CI: (0.52, 0.97)
mannwhitney(ord)
#> Mann-Whitney Estimand
#> 
#>      Estimator: tmle
#>       Estimate: 0.45
#>     Std. error: 0.01
#>         95% CI: (0.42, 0.47)
cdf(ord)
#> CDF Estimator: tmle
#> 
#> Arm-specific CDF: Pr(K <= k | A = a)
#> Treatment Arm
#>   k Estimate Std. error         95% CI Uniform 95% CI
#> 1 0    0.599      0.018 (0.56 to 0.63) (0.44 to 0.76)
#> 2 1    0.684      0.017 (0.65 to 0.72) (0.53 to 0.84)
#> 3 2    0.741      0.016 (0.71 to 0.77) (0.60 to 0.89)
#> 4 3    0.903      0.011 (0.88 to 0.92) (0.80 to 1.00)
#> 5 4    0.975      0.006 (0.96 to 0.99) (0.92 to 1.00)
#> 6 5    1.000          -              -              -
#> 
#> Control Arm
#>   k Estimate Std. error         95% CI Uniform 95% CI
#> 1 0    0.557      0.018 (0.52 to 0.59) (0.47 to 0.64)
#> 2 1    0.613      0.018 (0.58 to 0.65) (0.53 to 0.70)
#> 3 2    0.687      0.017 (0.65 to 0.72) (0.61 to 0.77)
#> 4 3    0.715      0.016 (0.68 to 0.75) (0.64 to 0.79)
#> 5 4    0.879      0.012 (0.86 to 0.90) (0.82 to 0.94)
#> 6 5    1.000          -              -              -
pmf(ord)
#> PMF Estimator: tmle
#> 
#> Arm-specific PMF: Pr(K = k | A = a)
#> Treatment Arm
#>   k Estimate Std. error         95% CI Uniform 95% CI
#> 1 0    0.599      0.018 (0.56 to 0.63) (0.44 to 0.76)
#> 2 1    0.084      0.010 (0.06 to 0.10) (0.00 to 0.18)
#> 3 2    0.057      0.009 (0.04 to 0.07) (0.00 to 0.13)
#> 4 3    0.162      0.013 (0.14 to 0.19) (0.04 to 0.28)
#> 5 4    0.072      0.009 (0.05 to 0.09) (0.00 to 0.16)
#> 6 5    0.025      0.006 (0.01 to 0.04) (0.00 to 0.08)
#> 
#> Control Arm
#>   k Estimate Std. error         95% CI Uniform 95% CI
#> 1 0    0.557      0.018 (0.52 to 0.59) (0.47 to 0.64)
#> 2 1    0.056      0.008 (0.04 to 0.07) (0.02 to 0.10)
#> 3 2    0.075      0.010 (0.06 to 0.09) (0.03 to 0.12)
#> 4 3    0.028      0.006 (0.02 to 0.04) (0.00 to 0.06)
#> 5 4    0.164      0.014 (0.14 to 0.19) (0.10 to 0.23)
#> 6 5    0.121      0.012 (0.10 to 0.14) (0.06 to 0.18)
```

# References

Díaz, I., E. Colantuoni, D. F. Hanley, and M. Rosenblum (2019). Improved
precision in the analysis of randomized trials with survival outcomes,
without assuming proportional hazards. Lifetime Data Analysis 25 (3),
439–468.

Díaz, I., Colantuoni, E. and Rosenblum, M. (2016), Enhanced precision in
the analysis of randomized trials with ordinal outcomes. Biom, 72:
422-431. <https://doi.org/10.1111/biom.12450>

Vermeulen, K., Thas, O., and Vansteelandt, S. (2015), Increasing the
power of the Mann‐Whitney test in randomized experiments through
flexible covariate adjustment. *Statist. Med.*, 34, pages 1012– 1030.
doi:
[10.1002/sim.6386](https://doi.org/10.1002/sim.6386 "Link to external resource: 10.1002/sim.6386")

Edward H. Kennedy (2019) Nonparametric Causal Effects Based on
Incremental Propensity Score Interventions, Journal of the American
Statistical Association, 114:526, 645-656, DOI:
10.1080/01621459.2017.1422737
