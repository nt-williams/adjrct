
Engine <- R6::R6Class(
  "Engine",
  public = list(
    engine = NULL,
    lrnrs_cens = NULL,
    lrnrs_hzrd = NULL,
    initialize = function(lrnrs_cens = NULL, lrnrs_hzrd = NULL) {
      self$lrnrs_cens <- lrnrs_cens
      self$lrnrs_hzrd <- lrnrs_hzrd
    },
    eval_engine = function(estimator) {
      if (estimator == "km") {
        self$engine <- "glm"
        if (!is.null(self$lrnrs_cens) || !is.null(self$lrnrs_hzrd)) {
          warning("sl3 learners supplied with Kaplan-Meier estimator. Using `glm()` instead.",
                  call. = FALSE)
        }
      } else {
        if (is.null(self$lrnrs_cens) || is.null(self$lrnrs_hzrd)) {
          self$engine <- "glm"
        } else {
          has_sl3 <- check_sl3()
          self$engine <- ifelse(has_sl3, "sl3", "glm")
          if (!has_sl3) warning("sl3 not found. Using `glm()`.", call. = FALSE)
        }
      }
      invisible(self)
    }
  )
)

nuis_dr <- function(self) {
  switch(self$engine$engine,
         glm = nuis_glm(self),
         sl3 = nuis_sl3(self))
}

nuis_sl3 <- function(self) {
  task_H <- new_sl3(self$at_risk_evnt(), "evnt", self$predictors(), "survrctId")
  task_R <- new_sl3(self$at_risk_cens(), "cens", self$predictors(), "survrctId")

  pred_H0 <- new_sl3(self$turn_off(), "evnt", self$predictors(), "survrctId")
  pred_H1 <- new_sl3(self$turn_on(), "evnt", self$predictors(), "survrctId")
  pred_R0 <- new_sl3(self$turn_off(), "cens", self$predictors(), "survrctId")
  pred_R1 <- new_sl3(self$turn_on(), "cens", self$predictors(), "survrctId")

  ensm_H <- new_ensemble(self$engine$lrnrs_hzrd)
  ensm_R <- new_ensemble(self$engine$lrnrs_cens)

  fit_H <- run_ensemble(ensm_H, task_H, envir = environment())
  fit_R <- run_ensemble(ensm_R, task_R, envir = environment())
  fit_A <- glm(self$formula_trt(), data = self$at_risk_trt(), family = "binomial")

  list(
    hzrd_off = predict_sl3(fit_H, pred_H0, envir = environment()),
    hzrd_on  = predict_sl3(fit_H, pred_H1, envir = environment()),
    cens_off = predict_sl3(fit_R, pred_R0, envir = environment()),
    cens_on  = predict_sl3(fit_R, pred_R1, envir = environment()),
    trt_off  = 1 - predict(fit_A, newdata = self$turn_on(), type = "response"),
    trt_on   = predict(fit_A, newdata = self$turn_on(), type = "response")
  )
}

nuis_glm <- function(self) {
  fit_H <- glm(self$formula_hzrd(), data = self$at_risk_evnt(), family = "binomial")
  fit_R <- glm(self$formula_cens(), data = self$at_risk_cens(), family = "binomial")
  fit_A <- glm(self$formula_trt(), data = self$at_risk_trt(), family = "binomial")

  list(
    hzrd_off = predict(fit_H, newdata = self$turn_off(), type = "response"),
    hzrd_on  = predict(fit_H, newdata = self$turn_on(), type = "response"),
    cens_off = predict(fit_R, newdata = self$turn_off(), type = "response"),
    cens_on  = predict(fit_R, newdata = self$turn_on(), type = "response"),
    trt_off  = 1 - predict(fit_A, newdata = self$turn_on(), type = "response"),
    trt_on   = predict(fit_A, newdata = self$turn_on(), type = "response")
  )
}
