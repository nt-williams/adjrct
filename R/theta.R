
rmst_eif <- function(meta, aux) {
  trt <- meta$get_var("trt")
  id <- meta$surv_data[["survrctId"]]
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
  se <- sqrt(var(eif) / meta$nobs)
  list(arm1            = theta1,
       arm0            = theta0,
       theta           = theta1 - theta0,
       eif             = eif,
       std.error       = se,
       theta.conf.low  = (theta1 - theta0) - qnorm(0.975)*se,
       theta.conf.high = (theta1 - theta0) + qnorm(0.975)*se)
}

survprob_eif <- function(estimator, meta, aux) {
  vals <- switch(estimator,
                 tmle = survprob_eif_tmle(meta, aux),
                 aipw = survprob_eif_aipw(meta, aux),
                 km = survprob_eif_tmle(meta, aux))
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
  id <- meta$surv_data[["survrctId"]]
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

survprob_eif_aipw <- function(meta, aux) {
  trt <- meta$get_var("trt")
  id <- meta$surv_data[["survrctId"]]
  dt1 <- sum_by_id(meta$risk_evnt * trt * aux$Z1 * (meta$surv_data[["evnt"]] - aux$LH), id)
  dt0 <- sum_by_id(meta$risk_evnt * (1 - trt) * aux$Z0 * (meta$surv_data[["evnt"]] - aux$LH), id)
  dw1 <- do_rbind(aux$S1, id)[meta$all_time == 1, aux$time]
  dw0 <- do_rbind(aux$S0, id)[meta$all_time == 1, aux$time]
  theta1 <- mean(dt1 + dw1) # TODO ASK IVAN WHY WE ARE ADDING THESE HERE!
  theta0 <- mean(dt0 + dw0)
  eif <- as.vector(dt1 - dt0 + dw1 - dw0)
  list(theta1 = theta1,
       theta0 = theta0,
       theta  = theta1 - theta0,
       eif    = eif)
}

simul_ci <- function(res, n) {
  cv <- simul::simul_crit(lapply(res, function(x) x$theta),
                          lapply(res, function(x) x$eif),
                          n)
  for (i in 1:length(res)) {
    res[[i]]$unif.conf.low <- res[[i]]$theta - cv*res[[i]]$std.error
    res[[i]]$unif.conf.high <- res[[i]]$theta + cv*res[[i]]$std.error
  }
  res$mbcv <- cv
  return(res)
}
