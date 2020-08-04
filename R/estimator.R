
rmst <- function(data, trt, status, baseline, time, id, coarsen = 1,
                 horizon, estimator = c("tmle", "aipw", "ipw", "km"),
                 learners_trt = NULL, learners_cens = NULL, learners_hazard = NULL) {

  # prepare data
  meta <- Meta$new(data, trt, status, baseline, id, time, horizon,
                   coarsen, estimator, learners_trt, learners_cens, learners_hazard)

  # generate nuisance values
  nuis <- estimate_nuisance(meta, estimator)

  # compute estimator
  compute_rmst(meta, nuis, estimator)
}
