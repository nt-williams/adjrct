
nuis_dr <- function(self) {
  task_H <- new_sl3(self$at_risk_evnt(), "evnt", self$predictors(), "survrctId")
  task_R <- new_sl3(self$at_risk_cens(), "cens", self$predictors(), "survrctId")
  task_A <- new_sl3(self$at_risk_trt(), self$trt, self$covar, "survrctId")

  pred_H0 <- new_sl3(self$turn_off(), "evnt", self$predictors(), "survrctId")
  pred_H1 <- new_sl3(self$turn_on(), "evnt", self$predictors(), "survrctId")
  pred_R0 <- new_sl3(self$turn_off(), "cens", self$predictors(), "survrctId")
  pred_R1 <- new_sl3(self$turn_on(), "cens", self$predictors(), "survrctId")
  pred_A1 <- new_sl3(self$turn_on(), self$trt, self$covar, "survrctId")

  ensm_H <- new_ensemble(self$lrnrs_hzrd)
  ensm_R <- new_ensemble(self$lrnrs_cens)
  ensm_A <- new_ensemble(self$lrnrs_trt)

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

nuis_ua <- function(self) {


}
