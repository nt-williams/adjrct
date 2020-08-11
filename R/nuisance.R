
estimate_nuisance <- function(meta, estimator) {
  switch(estimator,
         tmle = nuisance.dr(meta),
         aipw = nuisance.dr(meta),
         ipw  = nuisance.ua(meta),
         km   = nuisance.ua(meta))
}

nuisance.dr <- function(meta) {

  # create sl3 tasks
  task_H  <- sw(initiate_sl3_task(meta$data[meta$im == 1, ], "lm", c("m", meta$trt, meta$covar), "binomial", meta$id))
  task_R  <- sw(initiate_sl3_task(meta$data[meta$jm == 1, ], "rm", c("m", meta$trt, meta$covar), "binomial", meta$id))
  task_A  <- sw(initiate_sl3_task(meta$data[meta$m == 1, ], meta$trt, meta$covar, "binomial", meta$id))
  pred_H0 <- sw(initiate_sl3_task(turn_off(meta$data, meta$trt), "lm", c("m", meta$trt, meta$covar), "binomial", meta$id))
  pred_H1 <- sw(initiate_sl3_task(turn_on(meta$data, meta$trt), "lm", c("m", meta$trt, meta$covar), "binomial", meta$id))
  pred_R0 <- sw(initiate_sl3_task(turn_off(meta$data, meta$trt), "rm", c("m", meta$trt, meta$covar), "binomial", meta$id))
  pred_R1 <- sw(initiate_sl3_task(turn_on(meta$data, meta$trt), "rm", c("m", meta$trt, meta$covar), "binomial", meta$id))
  pred_A1 <- sw(initiate_sl3_task(turn_on(meta$data, meta$trt), meta$trt, meta$covar, "binomial", meta$id))

  # create learners
  ensm_H <- initiate_ensemble("binomial", meta$learners_hazard)
  ensm_R <- initiate_ensemble("binomial", meta$learners_cens)
  ensm_A <- initiate_ensemble("binomial", meta$learners_trt)

  # run sl3
  fit_H <- run_ensemble(ensm_H, task_H, envir = environment())
  fit_R <- run_ensemble(ensm_R, task_R, envir = environment())
  fit_A <- run_ensemble(ensm_A, task_A, envir = environment())

  # generate predictions
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

