
#' Fit Nuisance Parameters and Create Metadata
#'
#' @param data
#' @param trt
#' @param status
#' @param baseline
#' @param time
#' @param id
#' @param coarsen
#' @param estimator
#' @param learners_trt
#' @param learners_cens
#' @param learners_hazard
#'
#' @return
#' @export
#'
#' @examples
fit_nuisance <- function(data, trt, status, baseline, time, id,
                         coarsen = 1, estimator = c("tmle", "aipw", "km"),
                         lrnrs_trt = NULL, lrnrs_cens = NULL, lrnrs_hazard = NULL) {
  metadata <- Survival$new(data, trt, status, baseline, id,
                           time, coarsen, match.arg(estimator),
                           lrnrs_trt, lrnrs_cens, lrnrs_hazard)
  metadata$estimate_nuisance()
}


#' Estimate Restricted Mean Survival Time
#'
#' @param metadata
#' @param horizon
#'
#' @return
#' @export
#'
#' @examples
rmst <- function(metadata, horizon = NULL) {
  metadata$evaluate_horizon(horizon)
  out <- list(estimator = metadata$estimator,
              horizon   = metadata$horizon,
              estimates = compute_rmst(metadata))
  class(out) <- "rmst"
  out
}

#' Estimate Survival Probabilities
#'
#' @param metadata
#' @param horizon
#'
#' @return
#' @export
#'
#' @examples
survprob <- function(metadata, horizon = NULL) {
  metadata$evaluate_horizon(horizon)
  out <- list(estimator = estimator,
              horizon   = metadata$horizon,
              estimates = compute_survprob(metadata))
  class(out) <- "survprob"
  out
}
