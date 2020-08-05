
rmst_eif <- function(meta, trt, z1, z0, s1, s0, lh, id) {
  dt     <- sum_by_id(meta$im * (trt*z1 - (1 - trt)*z0) * (meta$data[["lm"]] - lh), id)
  dw1    <- rowSums(do.call('rbind', s1[id])[meta$m == 1, 1:(meta$tau - 1)])
  dw0    <- rowSums(do.call('rbind', s0[id])[meta$m == 1, 1:(meta$tau - 1)])
  theta1 <- 1 + mean(dw1)
  theta0 <- 1 + mean(dw0)
  theta  <- (1 + theta1) - (1 + theta0)
  eif    <- as.vector(dt + dw1 - dw0 - theta)
  se     <- sqrt(var(eif) / meta$n)

  list(rmst1     = theta1,
       rmst0     = theta0,
       theta     = theta,
       eif       = eif,
       std.error = se,
       theta.conf.low  = theta - qnorm(0.975)*se,
       theta.conf.high = theta + qnorm(0.975)*se)
}
