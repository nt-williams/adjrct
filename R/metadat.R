
Survival <- R6::R6Class(
  "Survival",
  public = list(
    data = NULL,
    estimator = NULL,
    trt = NULL,
    status = NULL,
    covar = NULL,
    time = NULL,
    horizon = list(),
    id = NULL,
    lrnrs_trt = NULL,
    lrnrs_cens = NULL,
    lrnrs_hazard = NULL,
    n = NULL,
    m = NULL,
    k = NULL,
    im = NULL,
    jm = NULL,
    nuisance = list(),
    initialize = function(data, trt, status, baseline, id, time,
                          coarsen, estimator, lrnrs_trt,
                          lrnrs_cens, lrnrs_hazard) {

      # checks
      # check_correct_trt() i.e., trt is 0 and 1
      # check_time_horizon() i.e., less than maximum time
      # check_status() i.e., status is 0 and 1

      tfd               <- prepare_data(data, time, status, coarsen)
      self$data         <- tfd$data
      self$trt          <- trt
      self$covar        <- baseline
      self$id           <- id
      self$time         <- time
      self$n            <- length(unique(tfd$data[[id]]))
      self$m            <- tfd$m
      self$k            <- max(tfd$m)
      self$im           <- tfd$im
      self$jm           <- tfd$jm
      self$estimator    <- estimator
      self$lrnrs_trt    <- check_sl3_usage("trt", estimator, lrnrs_trt)
      self$lrnrs_cens   <- check_sl3_usage("cens", estimator, lrnrs_cens)
      self$lrnrs_hazard <- check_sl3_usage("hazard", estimator, lrnrs_hazard)

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
        self$horizon <- 1:max(self$k)
      } else {
        self$horizon <- horizon
      }
      invisible(self)
    },
    print = function(...) {
      cat("rctSurv metadata")
    }
  )
)
