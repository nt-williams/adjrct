
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
> Proportional Hazards and Odds Assumptions

Nick Williams and Iván Díaz

-----

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
specifying the left-hand side of the formula.

We also specify the target variable of interest, an optional time
coarsener, and the estimator.

``` r
library(adjrct)

data(colon)
surv <- survrct(Surv(time, status) ~ trt + age + sex, target = "trt", data = colon, coarsen = 30)
#> survrct metadata
#> 
#> Surv(time, status) ~ trt + age + sex
#> 
#> ● Estimate RMST with `rmst()`
#> ● Estimate survival probability with `survprob()`
#> ● Inspect nuisance parameter models with `get_fits()`
#> 
#>          Estimator: tmle
#>    Target variable: trt
#>   Status Indicator: status
#>     Adjustment set: age and sex
#> Max coarsened time: 111
```

Using the metadata from the previous step we can now estimate the
restricted mean survival time and survival probability for a single or
multiple time horizons. If multiple times are evaluated, two confidence
bands are returned: 95% point-wise intervals as well as 95% uniform
confidence bands based on the multiplier-bootstrap from Kennedy (2019).

``` r
rmst(surv, 60)
#> RMST Estimator: tmle
#> 
#> Marginal RMST: E(min[T, 60] | A = a)
#> Treatment Arm
#>       Estimate: 45.33
#>     Std. error: 1.23
#>         95% CI: (42.91, 47.74)
#> Control Arm
#>       Estimate: 37.83
#>     Std. error: 1.39
#>         95% CI: (35.12, 40.55)
#> 
#> Treatment Effect: E(min[T, 60] | A = 1) - E(min[T, 60] | A = 0)
#> Additive effect
#>       Estimate: 7.49
#>     Std. error: 1.85
#>         95% CI: (3.86, 11.12)
survprob(surv, 60)
#> Survival Probability Estimator: tmle
#> 
#> Marginal Survival Probability: Pr(T > 60 | A = a)
#> Treatment Arm
#>       Estimate: 0.65
#>     Std. error: 0.03
#>         95% CI: (0.6, 0.7)
#> Control Arm
#>       Estimate: 0.51
#>     Std. error: 0.03
#>         95% CI: (0.45, 0.57)
#> 
#> Treatment Effect: Pr(T > 60 | A = 1) - Pr(T > 60 | A = 0)
#> Additive effect
#>       Estimate: 0.14
#>     Std. error: 0.04
#>         95% CI: (0.06, 0.22)
```

### Ordinal outcome

We can similarly create an Ordinal metadata object using the
`ordinalrct()` function.

``` r
data(mistie)

ord <- ordinalrct(Y ~ A + age, "A", mistie, "tmle", lasso = FALSE)
ord
#> ordinalrct metadata
#> 
#> Y ~ A + age
#> 
#> ● Estimate log odds ratio with `log_or()`
#> ● Estimate Mann-Whitney with `mannwhitney()`
#> ● Estimate with `cdf()`
#> ● Estimate with `pmf()`
#> ● Inspect nuisance parameter models with `get_fits()`
#> 
#>          Estimator: tmle
#>    Target variable: A
#>   Outcome variable: Y
#>     Adjustment set: age
```

The average log odds ratio, Mann-Whitney statistic, CDF, and PMF can
then be estimated using the metadata.

``` r
log_or(ord)
#> Log OR Estimator: tmle
#> 
#> Arm-specific log odds:
#> Treatment Arm
#>       Estimate: -0.24
#>     Std. error: 0.22
#>         95% CI: (-0.68, 0.19)
#> Control Arm
#>       Estimate: -0.5
#>     Std. error: 0.31
#>         95% CI: (-1.1, 0.11)
#> 
#> Average log odds ratio:
#>       Estimate: 0.25
#>     Std. error: 0.38
#>         95% CI: (-0.48, 0.99)
mannwhitney(ord)
#> Mann-Whitney Estimand
#> 
#>      Estimator: tmle
#>       Estimate: 0.46
#>     Std. error: 0.06
#>         95% CI: (0.34, 0.57)
cdf(ord)
#> CDF Estimator: tmle
#> 
#> Arm-specific CDF: Pr(K <= k | A = a)
#> Treatment Arm
#>   k Estimate Std. error         95% CI
#> 1 1    0.121      0.040 (0.04 to 0.20)
#> 2 2    0.364      0.059 (0.25 to 0.48)
#> 3 3    0.606      0.060 (0.49 to 0.72)
#> 4 4    0.757      0.053 (0.65 to 0.86)
#> 5 5    1.000          -              -
#> 
#> Control Arm
#>   k Estimate Std. error         95% CI
#> 1 1    0.108      0.051 (0.01 to 0.21)
#> 2 2    0.243      0.070 (0.11 to 0.38)
#> 3 3    0.566      0.082 (0.41 to 0.73)
#> 4 4    0.730      0.073 (0.59 to 0.87)
#> 5 5    1.000          -              -
pmf(ord)
#> PMF Estimator: tmle
#> 
#> Arm-specific PMF: Pr(K = k | A = a)
#> Treatment Arm
#>   k Estimate Std. error         95% CI
#> 1 1    0.121      0.040 (0.04 to 0.20)
#> 2 2    0.242      0.053 (0.14 to 0.35)
#> 3 3    0.242      0.053 (0.14 to 0.35)
#> 4 4    0.151      0.045 (0.06 to 0.24)
#> 5 5    0.243      0.053 (0.14 to 0.35)
#> 
#> Control Arm
#>   k Estimate Std. error         95% CI
#> 1 1    0.108      0.051 (0.01 to 0.21)
#> 2 2    0.135      0.056 (0.02 to 0.25)
#> 3 3    0.323      0.076 (0.17 to 0.47)
#> 4 4    0.163      0.061 (0.04 to 0.28)
#> 5 5    0.270      0.073 (0.13 to 0.41)
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
