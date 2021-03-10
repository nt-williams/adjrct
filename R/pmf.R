compute_pmf <- function(x) {
  switch(x$estimator,
         tmle = pmf_tmle(x))
}

calc_pmf <- function(x) {
  t(apply(x, 1, calc_pmf_k))
}

calc_pmf_k <- function(cdf) {
  K <- length(cdf) + 1
  pmf <- vector("numeric", K)
  for (i in K:2) {
    if (i == K) {
      pmf[i] <- 1 - cdf[i - 1]
    } else {
      pmf[i] <- cdf[i] - cdf[i - 1]
    }
  }
  pmf[1] <- cdf[1]
  return(pmf)
}

pmf_tmle <- function(meta) {
  CDF <- cdf_tmle(meta)
  PMF <- calc_pmf(CDF$dist)
  std.error <- t(sapply(CDF$eif, function(x) sqrt(diag(var(pmf_eif(x))) / meta$nobs)))
  ci <- list(theta1 = dist_ci(PMF[1, ], std.error[1, ]),
             theta0 = dist_ci(PMF[2, ], std.error[2, ]))
  list(dist = PMF,
       std.error = std.error,
       ci = ci,
       eif = lapply(CDF$eif, pmf_eif))
}

pmf_eif <- function(x) {
  out <- matrix(nrow = nrow(x), ncol = ncol(x) + 1)
  for (i in 2:(ncol(x) + 1)) {
    if (i == ncol(out)) {
      out[, i] <- 1 - x[, i - 1]
    } else {
      out[, i] <- x[, i] - x[, i - 1]
    }
  }
  out[, 1] <- x[, 1]
  return(out)
}
