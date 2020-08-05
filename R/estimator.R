
#' Title
#'
#' @param data
#' @param trt
#' @param status
#' @param baseline
#' @param time
#' @param id
#' @param horizon
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
rmst <- function(data, trt, status, baseline, time, id, horizon,
                 coarsen = 1, estimator = c("tmle", "aipw", "ipw", "km"),
                 learners_trt = NULL, learners_cens = NULL, learners_hazard = NULL) {

  # prepare data
  meta <- Meta$new(data, trt, status, baseline, id, time, horizon,
                   coarsen, estimator, learners_trt, learners_cens, learners_hazard)

  # generate nuisance values
  nuis <- estimate_nuisance(meta, estimator)

  # compute estimator
  psi <- compute_rmst(meta, nuis, estimator)

  # return estimates
  out <- list(estimator = estimator,
              horizon = horizon,
              estimates = psi)
  class(out) <- "rct_rmst"
  return(out)
}
