
Survival <- R6::R6Class(
  "Survival",
  public = list(
    data = NULL,
    surv_data = NULL,
    estimator = NULL,
    trt = NULL,
    status = NULL,
    covar = NULL,
    time = NULL,
    horizon = list(),
    id = NULL,
    lrnrs_trt = NULL,
    lrnrs_cens = NULL,
    lrnrs_hzrd = NULL,
    nobs = NULL,
    all_time = NULL,
    max_time = NULL,
    risk_evnt = NULL,
    risk_cens = NULL,
    nuisance = list(),
    initialize = function(data, trt, status, baseline, id, time,
                          estimator, lrnrs_trt, lrnrs_cens, lrnrs_hzrd) {

      # check_correct_trt() i.e., trt is 0 and 1
      # check_time_horizon() i.e., less than maximum time
      # check_status() i.e., status is 0 and 1

      self$data       <- data
      self$trt        <- trt
      self$covar      <- baseline
      self$status     <- status
      self$id         <- id
      self$time       <- time
      self$estimator  <- estimator
      self$lrnrs_trt  <- check_sl3_usage("trt", estimator, lrnrs_trt)
      self$lrnrs_cens <- check_sl3_usage("cens", estimator, lrnrs_cens)
      self$lrnrs_hzrd <- check_sl3_usage("hzrd", estimator, lrnrs_hzrd)

    },
    prepare_data = function(coarsen = 1) {
      if (coarsen > 1) self$data[[self$time]] <- self$data[[self$time]] %/% coarsen + 1

      nobs      <- nrow(self$data)
      max_time  <- max(self$data[[self$time]])
      all_time  <- rep(1:max_time, nobs)
      evnt      <- cens <- rep(NA, nobs*max_time)
      risk_evnt <- risk_cens <- 1*(all_time == 1)

      for (t in 1:max_time) {
        cens[all_time == t] <- (1 - self$data[[self$status]]) * (self$data[[self$time]] == t)
        evnt[all_time == t] <- self$data[[self$status]] * (self$data[[self$time]] == t)
        risk_evnt[all_time == t] <- (self$data[[self$time]] >= t)
        risk_cens[all_time == t] <- (self$data[[self$time]] > t) * self$data[[self$status]] +
          (self$data[[self$time]] >= t) * (1 - self$data[[self$status]])
      }

      self$surv_data <- data.frame(self$data[as.numeric(gl(nobs, max_time)), ], all_time, evnt, cens)
      self$nobs      <- nobs
      self$max_time  <- max_time
      self$all_time  <- all_time
      self$risk_evnt <- risk_evnt
      self$risk_cens <- risk_cens
      invisible(self)
    },
    estimate_nuisance = function(...) {
      self$nuisance <- switch(self$estimator,
                              tmle = nuisance.dr(self),
                              aipw = nuisance.dr(self),
                              km   = nuisance.ua(self))
      invisible(self)
    },
    evaluate_horizon = function(horizon = NULL) {
      if (is.null(horizon)) {
        self$horizon <- 1:self$max_time
      } else {
        self$horizon <- horizon
      }
      invisible(self)
    },
    turn_on = function() {
      out <- self$surv_data
      out[[self$trt]] <- rep(1, nrow(self$surv_data))
      return(out)
    },
    turn_off = function() {
      out <- self$surv_data
      out[[self$trt]] <- rep(0, nrow(self$surv_data))
      return(out)
    },
    get_var = function(var) {
      self$surv_data[[self[[var]]]]
    },
    time_indicator = function() {
      outer(self$all_time, 1:self$max_time, "<=")
    },
    print = function(...) {
      cat("rctSurv metadata \n")
      cli::cli_ul(c("Estimate RMST with `rmst()`",
                    "Estimate survival probability with `survprob()`",
                    "Inspect SuperLearner weights with `get_weights()`"))
    }
  )
)
