
nuisance.dr <- function(self) {
  task_H  <- sw(initiate_sl3_task(self$data[self$im == 1, ], "lm", c("m", self$trt, self$covar), "binomial", self$id))
  task_R  <- sw(initiate_sl3_task(self$data[self$jm == 1, ], "rm", c("m", self$trt, self$covar), "binomial", self$id))
  task_A  <- sw(initiate_sl3_task(self$data[self$m == 1, ], self$trt, self$covar, "binomial", self$id))
  pred_H0 <- sw(initiate_sl3_task(turn_off(self$data, self$trt), "lm", c("m", self$trt, self$covar), "binomial", self$id))
  pred_H1 <- sw(initiate_sl3_task(turn_on(self$data, self$trt), "lm", c("m", self$trt, self$covar), "binomial", self$id))
  pred_R0 <- sw(initiate_sl3_task(turn_off(self$data, self$trt), "rm", c("m", self$trt, self$covar), "binomial", self$id))
  pred_R1 <- sw(initiate_sl3_task(turn_on(self$data, self$trt), "rm", c("m", self$trt, self$covar), "binomial", self$id))
  pred_A1 <- sw(initiate_sl3_task(turn_on(self$data, self$trt), self$trt, self$covar, "binomial", self$id))

  ensm_H <- initiate_ensemble("binomial", self$lrnrs_hazard)
  ensm_R <- initiate_ensemble("binomial", self$lrnrs_cens)
  ensm_A <- initiate_ensemble("binomial", self$lrnrs_trt)

  fit_H <- run_ensemble(ensm_H, task_H, envir = environment())
  fit_R <- run_ensemble(ensm_R, task_R, envir = environment())
  fit_A <- run_ensemble(ensm_A, task_A, envir = environment())

  list(
    hazard_off = predict_sl3(fit_H, pred_H0, envir = environment()),
    hazard_on  = predict_sl3(fit_H, pred_H1, envir = environment()),
    cens_off   = predict_sl3(fit_R, pred_R0, envir = environment()),
    cens_on    = predict_sl3(fit_R, pred_R1, envir = environment()),
    treat_off  = 1 - predict_sl3(fit_A, pred_A1, envir = environment()),
    treat_on   = predict_sl3(fit_A, pred_A1, envir = environment())
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
      lm ~ 0 + offset(qlogis(lh)) + I(A * z1) + I((1 - A) * z0),
      family = binomial,
      subset = im == 1,
      data = data
    )
  ))
}

tilt_params.gamma <- function(data) {
  check_na_coef(coef(
    glm2::glm2(
      rm ~ 0 + offset(qlogis(gr)) + h,
      family = binomial,
      subset = jm == 1,
      data = data
    )
  ))
}

tilt_params.nu <- function(data) {
  check_na_coef(coef(
    glm2::glm2(
      A ~ 0 + offset(qlogis(a1)) + M,
      family = binomial,
      subset = m == 1,
      data = data
    )
  ))
}

