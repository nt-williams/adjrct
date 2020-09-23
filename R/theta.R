
rmst_eif <- function(meta, aux) {
  trt <- meta$get_var("trt")
  id <- meta$surv_data[["survrctId"]]

  DT1 <- sum_by_id(meta$risk_evnt*trt*aux$Z1*(meta$surv_data[["evnt"]] - aux$LH), id)
  DT0 <- sum_by_id(meta$risk_evnt*(1 - trt)*aux$Z0*(meta$surv_data[["evnt"]] - aux$LH), id)
  DW1 <- rowSums(matrix(aux$ST1[meta$all_time == 1, 1:(aux$time - 1)],
                        nrow = length(unique(id)),
                        ncol = aux$time - 1))
  DW0 <- rowSums(matrix(aux$ST0[meta$all_time == 1, 1:(aux$time - 1)],
                        nrow = length(unique(id)),
                        ncol = aux$time - 1))

  theta1 <- 1 + mean(DW1)
  theta0 <- 1 + mean(DW0)

  eif1 <- as.vector(DT1 + DW1)
  eif0 <- as.vector(DT0 + DW0)
  eif <- eif1 - eif0

  se1 <- sqrt(var(eif1) / meta$nobs)
  se0 <- sqrt(var(eif0) / meta$nobs)
  se <- sqrt(var(eif) / meta$nobs)

  list(arm1            = theta1,
       eif1            = eif1,
       arm1.std.error  = se1,
       arm1.conf.low   = theta1 - qnorm(0.975)*se1,
       arm1.conf.high  = theta1 + qnorm(0.975)*se1,
       arm0            = theta0,
       eif0            = eif0,
       arm0.std.error  = se0,
       arm0.conf.low   = theta0 - qnorm(0.975)*se0,
       arm0.conf.high  = theta0 + qnorm(0.975)*se0,
       theta           = theta1 - theta0,
       eif             = eif,
       std.error       = se,
       theta.conf.low  = (theta1 - theta0) - qnorm(0.975)*se,
       theta.conf.high = (theta1 - theta0) + qnorm(0.975)*se)
}

# survprob_eif <- function(estimator, meta, aux) {
#   vals <- switch(estimator,
#                  tmle = survprob_eif_tmle(meta, aux),
#                  aipw = survprob_eif_aipw(meta, aux),
#                  km = survprob_eif_tmle(meta, aux))
#   se <- sqrt(var(vals$eif) / meta$nobs)
#   list(arm1            = vals$theta1,
#        arm0            = vals$theta0,
#        theta           = vals$theta,
#        eif             = vals$eif,
#        std.error       = se,
#        theta.conf.low  = vals$theta - qnorm(0.975)*se,
#        theta.conf.high = vals$theta + qnorm(0.975)*se)
# }
#
# survprob_eif_tmle <- function(meta, aux) {
#   trt <- meta$get_var("trt")
#   id <- meta$surv_data[["survrctId"]]
#   dt <- sum_by_id(meta$risk_evnt * (trt*aux$Z1 - (1 - trt)*aux$Z0) * (meta$surv_data[["evnt"]] - aux$LH), id)
#   dw1 <- do_rbind(aux$S1, id)[meta$all_time == 1, aux$time]
#   dw0 <- do_rbind(aux$S0, id)[meta$all_time == 1, aux$time]
#   theta1 <- mean(dw1)
#   theta0 <- mean(dw0)
#   eif <- as.vector(dt + dw1 - dw0)
#   list(theta1 = theta1,
#        theta0 = theta0,
#        theta  = theta1 - theta0,
#        eif    = eif)
# }

survprob_eif <- function(meta, aux) {
  trt <- meta$get_var("trt")
  id <- meta$surv_data[["survrctId"]]

  DT1 <- sum_by_id(meta$risk_evnt*trt*aux$Z1*(meta$surv_data[["evnt"]] - aux$LH), id)
  DT0 <- sum_by_id(meta$risk_evnt*(1 - trt)*aux$Z0*(meta$surv_data[["evnt"]] - aux$LH), id)
  DW1 <- aux$ST1[meta$all_time == 1, aux$time]
  DW0 <- aux$ST0[meta$all_time == 1, aux$time]

  theta1 <- mean(DW1)
  theta0 <- mean(DW0)

  eif1 <- as.vector(DT1 + DW1)
  eif0 <- as.vector(DT0 + DW0)
  eif <- eif1 - eif0

  se1 <- sqrt(var(eif1) / meta$nobs)
  se0 <- sqrt(var(eif0) / meta$nobs)
  se <- sqrt(var(eif) / meta$nobs)

  list(arm1            = theta1,
       eif1            = eif1,
       arm1.std.error  = se1,
       arm1.conf.low   = theta1 - qnorm(0.975)*se1,
       arm1.conf.high  = theta1 + qnorm(0.975)*se1,
       arm0            = theta0,
       eif0            = eif0,
       arm0.std.error  = se0,
       arm0.conf.low   = theta0 - qnorm(0.975)*se0,
       arm0.conf.high  = theta0 + qnorm(0.975)*se0,
       theta           = theta1 - theta0,
       eif             = eif,
       std.error       = se,
       theta.conf.low  = (theta1 - theta0) - qnorm(0.975)*se,
       theta.conf.high = (theta1 - theta0) + qnorm(0.975)*se)
}

# survprob_eif_aipw <- function(meta, aux) {
#   trt <- meta$get_var("trt")
#   id <- meta$surv_data[["survrctId"]]
#   dt1 <- sum_by_id(meta$risk_evnt * trt * aux$Z1 * (meta$surv_data[["evnt"]] - aux$LH), id)
#   dt0 <- sum_by_id(meta$risk_evnt * (1 - trt) * aux$Z0 * (meta$surv_data[["evnt"]] - aux$LH), id)
#   dw1 <- do_rbind(aux$S1, id)[meta$all_time == 1, aux$time]
#   dw0 <- do_rbind(aux$S0, id)[meta$all_time == 1, aux$time]
#   theta1 <- mean(dt1 + dw1) # TODO ASK IVAN WHY WE ARE ADDING THESE HERE!
#   theta0 <- mean(dt0 + dw0)
#   eif <- as.vector(dt1 - dt0 + dw1 - dw0)
#   list(theta1 = theta1,
#        theta0 = theta0,
#        theta  = theta1 - theta0,
#        eif    = eif)
# }

simul_ci <- function(res, n) {
  cv_theta <- simul::simul_crit(lapply(res, function(x) x$theta),
                                lapply(res, function(x) x$eif), n)
  cv_arm1 <- simul::simul_crit(lapply(res, function(x) x$arm1),
                               lapply(res, function(x) x$eif1), n)
  cv_arm0 <- simul::simul_crit(lapply(res, function(x) x$arm0),
                               lapply(res, function(x) x$eif0), n)
  for (i in 1:length(res)) {
    res[[i]]$arm1.unif.low <- res[[i]]$arm1 - cv_arm1*res[[i]]$arm1.std.error
    res[[i]]$arm1.unif.high <- res[[i]]$arm1 + cv_arm1*res[[i]]$arm1.std.error
    res[[i]]$arm0.unif.low <- res[[i]]$arm0 - cv_arm0*res[[i]]$arm0.std.error
    res[[i]]$arm0.unif.high <- res[[i]]$arm0 + cv_arm0*res[[i]]$arm0.std.error
    res[[i]]$theta.unif.low <- res[[i]]$theta - cv_theta*res[[i]]$std.error
    res[[i]]$theta.unif.high <- res[[i]]$theta + cv_theta*res[[i]]$std.error
  }
  res$mbcv_theta <- cv_theta
  res$mbcv_treatment <- cv_arm1
  res$mbcv_control <- cv_arm0
  return(res)
}
