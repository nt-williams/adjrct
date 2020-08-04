
Meta <- R6::R6Class(
  "Meta",
  public = list(
    data = NULL,
    trt = NULL,
    status = NULL,
    covar = NULL,
    horizon = NULL,
    id = NULL,
    learners_trt = NULL,
    learners_cens = NULL,
    learners_hazard = NULL,
    n = NULL,
    m = NULL,
    im = NULL,
    jm = NULL,
    initialize = function(data, trt, status, baseline, id, time, horizon,
                          coarsen, estimator, learners_trt, learners_cens, learners_hazard) {

      # checks

      # data prep
      tfd       <- prepare_data(data, time, status, coarsen)
      self$data <- tfd$data

      # necessary parameters
      self$trt     <- trt
      self$covar   <- baseline
      self$id      <- id
      self$horizon <- horizon
      self$n       <- tfd$n
      self$m       <- tfd$m
      self$im      <- tfd$im
      self$jm      <- tfd$jm

      # initiate sl3 if too be used
      self$learners_trt    <- check_sl3_usage("trt", estimator, learners_trt)
      self$learners_cens   <- check_sl3_usage("cens", estimator, learners_cens)
      self$learners_hazard <- check_sl3_usage("hazard", estimator, learners_hazard)

    }
  )
)

prepare_data <- function(data, time, status, coarsen) {

  if (coarsen > 1) data[[time]] <- data[[time]] %/% coarsen + 1

  n <- nrow(data)
  k <- max(data[[time]])
  m <- rep(1:k, n)
  lm <- rm <- rep(NA, n*k)
  im <- jm <- 1*(m == 1)

  for (t in 1:k) {
    rm[m == t] <- (1 - data[[status]]) * (data[[time]] == t)
    lm[m == t] <- data[[status]] * (data[[time]] == t)
    im[m == t] <- (data[[time]] >= t)
    jm[m == t] <- (data[[time]] > t) * data[[status]] +
      (data[[time]] >= t) * (1 - data[[status]])
  }

  list(data = data.frame(data[as.numeric(gl(n, k)), ], lm, rm),
       n = n,
       m = m,
       im = im,
       jm = jm)

}

turn_on <- function(data, trt) {
  data[[trt]] <- rep(1, nrow(data))
  return(data)
}

turn_off <- function(data, trt) {
  data[[trt]] <- rep(0, nrow(data))
  return(data)
}
