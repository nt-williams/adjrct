compute_mw <- function(x) {
  switch(x$estimator,
         tmle = mw_tmle(x))
}

calc_mw <- function(CDF, PMF) {
  F_0 <- c(0, CDF$dist["theta0", ])
  f_0 <- PMF$dist["theta0", ]
  f_1 <- PMF$dist["theta1", ]
  sum((F_0 + 0.5*f_0) * f_1)
}

mw_tmle <- function(meta) {
  CDF <- cdf_tmle(meta)
  PMF <- pmf_tmle(meta)
  mw <- calc_mw(CDF, PMF)
  eif <- cbind(CDF$eif$theta0, PMF$eif$theta1, PMF$eif$theta0)
  grad <- mw_gradient(CDF, PMF)
  std.error <- sqrt(t(grad) %*% cov(eif) %*% grad / meta$nobs)
  list(mann.whitney = mw,
       std.error = std.error,
       ci = mw + qnorm(c(0.05 / 2, 1 - 0.05 / 2)) * c(std.error),
       eif = eif)
}

# based on https://github.com/benkeser/drord/blob/master/R/mannwhitney_fn.R
mw_gradient <- function(CDF, PMF) {
  F_0 <- c(0, CDF$dist["theta0", ])
  f_0 <- PMF$dist["theta0", ]
  f_1 <- PMF$dist["theta1", ]
  c(f_1[-1], F_0 + 0.5*f_0, 0.5*f_1)
}
