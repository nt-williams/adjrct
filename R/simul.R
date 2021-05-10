simul_ci <- function(res, n) {
  if (length(res) == 1) {
    res$mbcv_control <- res$mbcv_treatment <- res$mbcv_theta <- qnorm(0.975)
    return(res)
  }
  cv <- simul::simul(lapply(res, function(x) x$theta), lapply(res, function(x) x$eif), n, 1e4)
  cv1 <- simul::simul(lapply(res, function(x) x$arm1), lapply(res, function(x) x$eif1), n, 1e4)
  cv0 <- simul::simul(lapply(res, function(x) x$arm0), lapply(res, function(x) x$eif0), n, 1e4)
  for (i in 1:length(res)) {
    res[[i]]$arm1.unif.low <- res[[i]]$arm1 - cv1*res[[i]]$arm1.std.error
    res[[i]]$arm1.unif.high <- res[[i]]$arm1 + cv1*res[[i]]$arm1.std.error
    res[[i]]$arm0.unif.low <- res[[i]]$arm0 - cv0*res[[i]]$arm0.std.error
    res[[i]]$arm0.unif.high <- res[[i]]$arm0 + cv0*res[[i]]$arm0.std.error
    res[[i]]$theta.unif.low <- res[[i]]$theta - cv*res[[i]]$std.error
    res[[i]]$theta.unif.high <- res[[i]]$theta + cv*res[[i]]$std.error
  }
  res$mbcv_theta <- cv
  res$mbcv_treatment <- cv1
  res$mbcv_control <- cv0
  return(res)
}
