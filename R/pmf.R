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
  eif <- lapply(CDF$eif, pmf_eif)
  std.error <- t(sapply(eif, function(x) sqrt(diag(var(x)) / meta$nobs)))

  ci <- list(theta1 = dist_ci(PMF[1, ], std.error[1, ]),
             theta0 = dist_ci(PMF[2, ], std.error[2, ]),
             unif1 = dist_ci(PMF[1, ], std.error[1, ], CDF$mbcv["theta1"]),
             unif0 = dist_ci(PMF[2, ], std.error[2, ], CDF$mbcv["theta0"]))

  list(dist = PMF,
       std.error = std.error,
       ci = ci,
       eif = eif,
       mbcv = CDF$mbcv)
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

dist_data <- function(obj, type = c("cdf", "pmf")) {
  if (match.arg(type) == "cdf") {
    out <- data.frame(
      arm = rep(c("trt", "control"), each = length(obj$levels)),
      k = rep(ordered(obj$levels, levels = obj$levels), 2),
      estimate = c(obj$estimates$dist[1, ], 1, obj$estimates$dist[2, ], 1),
      std.error = c(obj$estimates$std.error[1, ], NA_real_, obj$estimates$std.error[2, ], NA_real_),
      conf.low = c(obj$estimates$ci$theta1[, 1], NA_real_, obj$estimates$ci$theta0[, 1], NA_real_),
      conf.high = c(obj$estimates$ci$theta1[, 2], NA_real_, obj$estimates$ci$theta0[, 2], NA_real_),
      conf.unif.low = c(obj$estimates$ci$unif1[, 1], NA_real_, obj$estimates$ci$unif0[, 1], NA_real_),
      conf.unif.high = c(obj$estimates$ci$unif1[, 2], NA_real_, obj$estimates$ci$unif0[, 2], NA_real_)
    )
    return(out)
  }

  data.frame(
    arm = rep(c("trt", "control"), each = length(obj$levels)),
    k = rep(ordered(obj$levels, levels = obj$levels), 2),
    estimate = c(obj$estimates$dist[1, ], obj$estimates$dist[2, ]),
    std.error = c(obj$estimates$std.error[1, ], obj$estimates$std.error[2, ]),
    conf.low = c(obj$estimates$ci$theta1[, 1], obj$estimates$ci$theta0[, 1]),
    conf.high = c(obj$estimates$ci$theta1[, 2], obj$estimates$ci$theta0[, 2]),
    conf.unif.low = c(obj$estimates$ci$unif1[, 1], obj$estimates$ci$unif0[, 1]),
    conf.unif.high = c(obj$estimates$ci$unif1[, 2], obj$estimates$ci$unif0[, 2])
  )
}
