
compute_simulband <- function(res, n) {
  cv <- compute_multboot(vapply(res, function(x) x$theta, FUN.VALUE = 1),
                         lapply(res, function(x) x$eif),
                         n,
                         10000)

  for (i in 1:length(res)) {
    res[[i]]$unif.conf.low <- res[[i]]$theta - cv*res[[i]]$std.error
    res[[i]]$unif.conf.high <- res[[i]]$theta + cv*res[[i]]$std.error
  }

  res$mbcv <- cv
  res
}


compute_multboot <- function(psi, eif, n, reps) {
  psi <- matrix(rep(psi, n), nrow = n, byrow = TRUE)
  sig <- matrix(rep(vapply(eif, function(e) sqrt(var(e)), FUN.VALUE = 1), n), nrow = n, byrow = TRUE)
  eif <- matrix(unlist(eif), nrow = n, byrow = FALSE)
  mbs <- matrix(2 * rbinom(n * reps, 1, 0.5) - 1, nrow = n, ncol = reps)
  sup <- sapply(1:reps, function(i) {
    max(abs(apply(mbs[, i] * ((eif - psi) / sig), 2, sum) / sqrt(n)))
  })
  as.vector(quantile(sup, 0.95))
}
