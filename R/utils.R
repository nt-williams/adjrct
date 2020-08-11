
Meta <- R6::R6Class(
  "Meta",
  public = list(
    data = NULL,
    trt = NULL,
    status = NULL,
    covar = NULL,
    time = NULL,
    tau = NULL,
    id = NULL,
    learners_trt = NULL,
    learners_cens = NULL,
    learners_hazard = NULL,
    n = NULL,
    m = NULL,
    k = NULL,
    im = NULL,
    jm = NULL,
    initialize = function(data, trt, status, baseline, id, time, horizon,
                          coarsen, estimator, learners_trt, learners_cens, learners_hazard) {

      # checks
      # check_correct_trt() i.e., trt is 0 and 1
      # check_time_horizon() i.e., less than maximum time
      # check_status() i.e., status is 0 and 1

      # data prep
      tfd       <- prepare_data(data, time, status, coarsen)
      self$data <- tfd$data

      # necessary parameters
      self$trt     <- trt
      self$covar   <- baseline
      self$id      <- id
      self$time    <- time
      self$tau     <- evaluate_horizon(max(tfd$m), horizon)
      self$n       <- length(unique(tfd$data[[id]]))
      self$m       <- tfd$m
      self$k       <- max(tfd$m)
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

  list(data = data.frame(data[as.numeric(gl(n, k)), ], m, lm, rm),
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

prodlag <- function(x) {
  cumprod(c(1, x[-length(x)]))
}

cumprod_by_id <- function(x, id) {
  tapply(x, id, cumprod, simplify = FALSE)
}

prodlag_by_id <- function(x, id) {
  tapply(x, id, prodlag, simplify = FALSE)
}

access_meta_var <- function(meta, var) {
  meta$data[[meta[[var]]]]
}

compute_z <- function(ind, time, id, S, A, R) {
  -rowSums(matrix((ind * do.call("rbind", S[id]))[, 1:(time - 1)],
                  nrow = length(id),
                  ncol = time - 1)) /
    bound(unlist(S) * A[id] * unlist(R))
}

compute_z_survprob <- function(ind, id, time, S, A, R) {
  -(ind * do.call("rbind", S[id]))[, time] /
    bound(unlist(S) * A[id] * unlist(R))
}

compute_h_survprob <- function(id, time, S, A, R) {
  -do.call("rbind", S[id])[, time] /
    bound(unlist(S) * A[id] * unlist(R))
}

sum_by_id <- function(x, id) {
  tapply(x, id, sum)
}

compute_m <- function(S1, S0, A1, A0, m, time, id) {
  (sum_by_id(unlist(S1) / bound(A1[id]) * (m <= time - 1), id) +
    sum_by_id(unlist(S0) / bound(A0[id]) * (m <= time - 1), id))[id]
}

compute_m_survprob <- function(S1, S0, A1, A0, m, time, id) {
  (sum_by_id(unlist(S1) / bound(A1[id]) * (m <= time), id) +
     sum_by_id(unlist(S0) / bound(A0[id]) * (m <= time), id))[id]
}

bound <- function(x, lower = 0.001){
  x[x < lower] <- lower
  return(as.numeric(x))
}

bound01 <- function(x, bound = 1e-10){
  x[x < bound] <- bound
  x[x > 1 - bound] <- 1 - bound
  return(as.numeric(x))
}

sw <- function(x) {
  suppressWarnings(x)
}

evaluate_horizon <- function(time, horizon) {
  if (is.null(horizon)) {
    1:max(time)
  } else {
    horizon
  }
}

get_workers <- function() {
  workers <- formals(plan("next"))$workers
  if (is.null(workers)) {
    workers <- 1
  }
  return(workers)
}

check_future_job <- function() {
  if (get_workers() == 1) {
    FALSE
  } else {
    TRUE
  }
}
