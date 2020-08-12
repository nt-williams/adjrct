
nuisance.dr <- function(self) {
  task_H  <- sw(initiate_sl3_task(self$surv_data[self$risk_evnt == 1, ], "evnt", c("all_time", self$trt, self$covar), "binomial", self$id))
  task_R  <- sw(initiate_sl3_task(self$surv_data[self$risk_cens == 1, ], "cens", c("all_time", self$trt, self$covar), "binomial", self$id))
  task_A  <- sw(initiate_sl3_task(self$surv_data[self$all_time == 1, ], self$trt, self$covar, "binomial", self$id))
  pred_H0 <- sw(initiate_sl3_task(self$turn_off(), "evnt", c("all_time", self$trt, self$covar), "binomial", self$id))
  pred_H1 <- sw(initiate_sl3_task(self$turn_on(), "evnt", c("all_time", self$trt, self$covar), "binomial", self$id))
  pred_R0 <- sw(initiate_sl3_task(self$turn_off(), "cens", c("all_time", self$trt, self$covar), "binomial", self$id))
  pred_R1 <- sw(initiate_sl3_task(self$turn_on(), "cens", c("all_time", self$trt, self$covar), "binomial", self$id))
  pred_A1 <- sw(initiate_sl3_task(self$turn_on(), self$trt, self$covar, "binomial", self$id))

  ensm_H <- initiate_ensemble("binomial", self$lrnrs_hzrd)
  ensm_R <- initiate_ensemble("binomial", self$lrnrs_cens)
  ensm_A <- initiate_ensemble("binomial", self$lrnrs_trt)

  fit_H <- run_ensemble(ensm_H, task_H, envir = environment())
  fit_R <- run_ensemble(ensm_R, task_R, envir = environment())
  fit_A <- run_ensemble(ensm_A, task_A, envir = environment())

  list(
    hzrd_off = predict_sl3(fit_H, pred_H0, envir = environment()),
    hzrd_on  = predict_sl3(fit_H, pred_H1, envir = environment()),
    cens_off = predict_sl3(fit_R, pred_R0, envir = environment()),
    cens_on  = predict_sl3(fit_R, pred_R1, envir = environment()),
    trt_off  = 1 - predict_sl3(fit_A, pred_A1, envir = environment()),
    trt_on   = predict_sl3(fit_A, pred_A1, envir = environment())
  )
}

nuisance.ua <- function(meta) {


}

tilt_params <- function(data, nuisance) {
  switch(nuisance,
         epsilon = tilt_params.epsilon(data),
         gamma   = tilt_params.gamma(data),
         nu      = tilt_params.nu(data))
}

tilt_params.epsilon <- function(data) {
  check_na_coef(coef(
    glm2::glm2(
      evnt ~ 0 + offset(qlogis(offset)) + I(trt * z1) + I((1 - trt) * z0),
      family = binomial,
      subset = risk_evnt == 1,
      data = data
    )
  ))
}

tilt_params.gamma <- function(data) {
  check_na_coef(coef(
    glm2::glm2(
      cens ~ 0 + offset(qlogis(offset)) + h,
      family = binomial,
      subset = risk_cens == 1,
      data = data
    )
  ))
}

tilt_params.nu <- function(data) {
  check_na_coef(coef(
    glm2::glm2(
      trt ~ 0 + offset(qlogis(offset)) + M,
      family = binomial,
      subset = time == 1,
      data = data
    )
  ))
}

