
#' Create RCT Survival Metadata
#'
#' @param data
#' @param trt
#' @param status
#' @param baseline
#' @param time
#' @param id
#' @param coarsen
#' @param estimator
#' @param lrnrs_trt
#' @param lrnrs_cens
#' @param lrnrs_hzrd
#'
#' @return
#' @export
#'
#' @examples
metadata <- function(data, trt, status, baseline, time, id,
                     coarsen = 1, estimator = c("tmle", "aipw", "km"),
                     lrnrs_trt = NULL, lrnrs_cens = NULL, lrnrs_hzrd = NULL) {
  Survival$
    new(data, trt, status, baseline, id, time, match.arg(estimator),
        lrnrs_trt, lrnrs_cens, lrnrs_hzrd)$
    prepare_data(coarsen)$
    fit_nuis()
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
