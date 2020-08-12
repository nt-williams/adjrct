
rmst_eif <- function(meta, estimator, trt, tau, z1, z0, s1, s0, lh, id) {

  vals <- switch(estimator,
                 tmle = rmst_eif.tmle(meta, trt, tau, z1, z0, s1, s0, lh, id),
                 aipw = rmst_eif.aipw(meta, trt, tau, z1, z0, s1, s0, lh, id),
                 ipw  = rmst_eif.ipw(meta, trt, z1, z0, id))

  se <- sqrt(var(vals$eif) / meta$nobs)

  list(arm1      = vals$theta1,
       arm0      = vals$theta0,
       theta     = vals$theta,
       eif       = vals$eif,
       std.error = se,
       theta.conf.low  = vals$theta - qnorm(0.975)*se,
       theta.conf.high = vals$theta + qnorm(0.975)*se)

}

rmst_eif.tmle <- function(meta, trt, tau, z1, z0, s1, s0, lh, id) {
  dt     <- sum_by_id(meta$risk_evnt * (trt*z1 - (1 - trt)*z0) * (meta$surv_data[["evnt"]] - lh), id)
  dw1    <- rowSums(matrix(do.call('rbind', s1[id])[meta$all_time == 1, 1:(tau - 1)],
                           nrow = length(unique(id)),
                           ncol = tau - 1))
  dw0    <- rowSums(matrix(do.call('rbind', s0[id])[meta$all_time == 1, 1:(tau - 1)],
                           nrow = length(unique(id)),
                           ncol = tau - 1))
  theta1 <- 1 + mean(dw1)
  theta0 <- 1 + mean(dw0)
  eif    <- as.vector(dt + dw1 - dw0)

  list(theta1 = theta1,
       theta0 = theta0,
       theta  = theta1 - theta0,
       eif    = eif)
}

rmst_eif.aipw <- function(meta, trt, z1, z0, s1, s0, lh, id) {
  dt1    <- sum_by_id(meta$im * trt*z1 * (meta$data[["lm"]] - lh), id)
  dt0    <- sum_by_id(meta$im * (1 - trt)*z0 * (meta$data[["lm"]] - lh), id)
  dw1    <- rowSums(do.call('rbind', s1[id])[meta$m == 1, 1:(meta$tau - 1)])
  dw0    <- rowSums(do.call('rbind', s0[id])[meta$m == 1, 1:(meta$tau - 1)])
  theta1 <- 1 + mean(dt1 + dw1)
  theta0 <- 1 + mean(dt0 + dw0)
  eif    <- as.vector(dt1 - dt0 + dw1 - dw0)

  list(theta1 = theta1,
       theta0 = theta0,
       theta  = theta1 - theta0,
       eif    = eif)
}

rmst_eif.ipw <- function(meta, trt, z1, z0, id) {
  dt1    <- sum_by_id(meta$im * trt * z1 * meta$data[["lm"]], id)
  dt0    <- sum_by_id(meta$im * (1 - trt)*z0 * meta$data[["lm"]], id)
  theta1 <- meta$tau + mean(dt1)
  theta0 <- meta$tau + mean(dt0)

  list(theta1 = theta1,
       theta0 = theta0,
       theta  = theta1 - theta0,
       eif    = NULL)
}

survprob_eif <- function(meta, estimator, trt, tau, z1, z0, s1, s0, lh, id) {

  vals <- switch(estimator,
                 tmle = survprob_eif.tmle(meta, trt, tau, z1, z0, s1, s0, lh, id))

  se <- sqrt(var(vals$eif) / meta$n)

  list(arm1            = vals$theta1,
       arm0            = vals$theta0,
       theta           = vals$theta,
       eif             = vals$eif,
       std.error       = se,
       theta.conf.low  = vals$theta - qnorm(0.975)*se,
       theta.conf.high = vals$theta + qnorm(0.975)*se)

}

survprob_eif.tmle <- function(meta, trt, tau, z1, z0, s1, s0, lh, id) {
  dt     <- sum_by_id(meta$im * (trt*z1 - (1 - trt)*z0) * (meta$data[["lm"]] - lh), id)
  dw1    <- do.call('rbind', s1[id])[meta$m == 1, tau]
  dw0    <- do.call('rbind', s0[id])[meta$m == 1, tau]
  theta1 <- mean(dw1)
  theta0 <- mean(dw0)
  eif    <- as.vector(dt + dw1 - dw0)

  list(theta1 = theta1,
       theta0 = theta0,
       theta  = theta1 - theta0,
       eif    = eif)
}
