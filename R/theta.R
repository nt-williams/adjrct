
rmst_eif <- function(estimator, meta, aux) {
  vals <- switch(estimator,
                 tmle = rmst_eif_dr(meta, aux),
                 aipw = rmst_eif_dr(meta, aux))
  se <- sqrt(var(vals$eif) / meta$nobs)
  list(arm1            = vals$theta1,
       arm0            = vals$theta0,
       theta           = vals$theta,
       eif             = vals$eif,
       std.error       = se,
       theta.conf.low  = vals$theta - qnorm(0.975)*se,
       theta.conf.high = vals$theta + qnorm(0.975)*se)
}

rmst_eif_dr <- function(meta, aux) {
  trt <- meta$get_var("trt")
  id <- meta$surv_data[["rctSurvId"]]
  dt <- sum_by_id(meta$risk_evnt * (trt*aux$Z1 - (1 - trt)*aux$Z0) * (meta$surv_data[["evnt"]] - aux$LH), id)
  dw1 <- rowSums(matrix(do_rbind(aux$S1, id)[meta$all_time == 1, 1:(aux$time - 1)],
                        nrow = length(unique(id)),
                        ncol = aux$time - 1))
  dw0 <- rowSums(matrix(do_rbind(aux$S0, id)[meta$all_time == 1, 1:(aux$time - 1)],
                           nrow = length(unique(id)),
                           ncol = aux$time - 1))
  theta1 <- 1 + mean(dw1)
  theta0 <- 1 + mean(dw0)
  eif <- as.vector(dt + dw1 - dw0)
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

survprob_eif <- function(estimator, meta, aux) {
  vals <- switch(estimator,
                 tmle = survprob_eif_tmle(meta, aux))
  se <- sqrt(var(vals$eif) / meta$nobs)
  list(arm1            = vals$theta1,
       arm0            = vals$theta0,
       theta           = vals$theta,
       eif             = vals$eif,
       std.error       = se,
       theta.conf.low  = vals$theta - qnorm(0.975)*se,
       theta.conf.high = vals$theta + qnorm(0.975)*se)
}

survprob_eif_tmle <- function(meta, aux) {
  trt <- meta$get_var("trt")
  id <- meta$surv_data[["rctSurvId"]]
  dt <- sum_by_id(meta$risk_evnt * (trt*aux$Z1 - (1 - trt)*aux$Z0) * (meta$surv_data[["evnt"]] - aux$LH), id)
  dw1 <- do_rbind(aux$S1, id)[meta$all_time == 1, aux$time]
  dw0 <- do_rbind(aux$S0, id)[meta$all_time == 1, aux$time]
  theta1 <- mean(dw1)
  theta0 <- mean(dw0)
  eif <- as.vector(dt + dw1 - dw0)
  list(theta1 = theta1,
       theta0 = theta0,
       theta  = theta1 - theta0,
       eif    = eif)
}
