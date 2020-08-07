
check_name_conflicts <- function(data) {


}

check_sl3_usage <- function(param, estimator, learners) {
  switch(param,
         cens   = check_sl3_usage.cens(estimator, learners),
         trt    = check_sl3_usage.trt(estimator, learners),
         hazard = check_sl3_usage.hazard(estimator, learners))
}

check_sl3_usage.trt <- function(estimator, learners) {
  if (estimator %in% c("km")) {
    asymp_warning("learners_trt", estimator)
    return(NULL)
  } else {
    learners
  }
}

check_sl3_usage.cens <- function(estimator, learners) {
  if (estimator %in% c("km")) {
    asymp_warning("learners_cens", estimator)
    return(NULL)
  } else {
    learners
  }
}

check_sl3_usage.hazard <- function(estimator, learners) {
  if (estimator %in% c("km")) {
    asymp_warning("learners_hazard", estimator)
    return(NULL)
  } else {
    learners
  }
}

asymp_warning <- function(param, estimator) {
  warning(param, " provided with estimator = ", estimator, ". Defaulting back to GLM.", call. = FALSE)
}

check_na_coef <- function(coefs) {
  coefs[is.na(coefs)] <- 0
  return(coefs)
}
