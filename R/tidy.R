#' @importFrom generics tidy
#' @export
generics::tidy

#' Tidy a(n) rmst object
#'
#' @param x A `rmst` object.
#' @param ... Unused, included for generic consistency only.
#' @return A tidy [tibble::tibble()] summarizing restricted mean survival time estimates.
#'
#' @examples
#' surv <- survrct(Surv(days, event) ~ A + age + sex + dyspnea + bmi,
#'                 A ~ 1, data = c19.tte)
#' tidy(rmst(surv, 14))
#'
#' @export
tidy.rmst <- function(x, ...) {
  out <- x$estimates[-which(names(x$estimates) %in% c("mbcv_theta", "mbcv_treatment", "mbcv_control"))]
  out <- do.call(
    "rbind",
    lapply(out, function(x) as.data.frame(x[-which(names(x) %in% c("eif", "eif1", "eif0"))]))
  )
  out <- cbind(horizon = x$horizon, out)
  names(out) <- gsub("arm0", "control.rmst", gsub("arm1", "trt.rmst", names(out)))
  class(out) <- c("tbl_df", "tbl", "data.frame")
  return(out)
}

#' Tidy a(n) survprob object
#'
#' @param x A `survprob` object.
#' @param ... Unused, included for generic consistency only.
#' @return A tidy [tibble::tibble()] summarizing survival probability estimates.
#'
#' @examples
#' surv <- survrct(Surv(days, event) ~ A + age + sex + dyspnea + bmi,
#'                 A ~ 1, data = c19.tte)
#' tidy(survprob(surv, 14))
#'
#' @export
tidy.survprob <- function(x, ...) {
  out <- tidy.rmst(x)
  names(out) <- gsub("rmst", "survprob", names(out))
  return(out)
}

#' Tidy a(n) lor object
#'
#' @param x A `lor` object.
#' @param ... Unused, included for generic consistency only.
#' @return A tidy [tibble::tibble()] summarizing log odds estimates.
#'
#' @examples
#' rct <- ordinalrct(state_ordinal ~ A + age, A ~ 1, data = c19.ordinal)
#' tidy(log_or(rct))
#'
#' @export
tidy.lor <- function(x, ...) {
  tmp <- x$estimates[-which(names(x$estimates) %in% c("eif", "lo1.ci", "lo0.ci", "ci"))]
  cis <- lapply(x$estimates[c("lo1.ci", "lo0.ci", "ci")], function(x) {
    names(x) <- c("conf.low", "conf.high")
    as.list(x)
  })
  out <- as.data.frame(c(tmp, as.list(unlist(cis))))
  names(out) <- gsub("ci.", "", gsub("lo0", "control.logodds", gsub("lo1", "trt.logodds", names(out))))
  class(out) <- c("tbl_df", "tbl", "data.frame")
  return(out)
}

#' Tidy a(n) cdf object
#'
#' @param x A `cdf` object.
#' @param ... Unused, included for generic consistency only.
#' @return A tidy [tibble::tibble()] summarizing CDF estimates.
#'
#' @examples
#' rct <- ordinalrct(state_ordinal ~ A + age, A ~ 1, data = c19.ordinal)
#' tidy(cdf(rct))
#'
#' @export
tidy.cdf <- function(x, ...) {
  dist <- t(x$estimates$dist)
  std.error <- x$estimates$std.error
  rownames(std.error) <- c("theta1.std.error", "theta0.std.error")
  ci <- x$estimates$ci
  ci <- mapply(function(x, n) {
    colnames(x) <- paste(n, c("conf.low", "conf.high"), sep = ".")
    x
  }, ci, names(ci), SIMPLIFY = FALSE)
  out <- cbind(k = x$levels[-length(x$levels)], as.data.frame(dist),
               as.data.frame(t(std.error)), as.data.frame(do.call("cbind", ci)))
  names(out) <- gsub("theta0", "control.cdf", gsub("theta1", "trt.cdf", names(out)))
  class(out) <- c("tbl_df", "tbl", "data.frame")
  return(out)
}

#' Tidy a(n) pmf object
#'
#' @param x A `pmf` object.
#' @param ... Unused, included for generic consistency only.
#' @return A tidy [tibble::tibble()] summarizing PMF estimates.
#'
#' @examples
#' rct <- ordinalrct(state_ordinal ~ A + age, A ~ 1, data = c19.ordinal)
#' tidy(pmf(rct))
#'
#' @export
tidy.pmf <- function(x, ...) {
  dist <- t(x$estimates$dist)
  std.error <- x$estimates$std.error
  rownames(std.error) <- c("theta1.std.error", "theta0.std.error")
  ci <- x$estimates$ci
  ci <- mapply(function(x, n) {
    colnames(x) <- paste(n, c("conf.low", "conf.high"), sep = ".")
    x
  }, ci, names(ci), SIMPLIFY = FALSE)
  out <- cbind(k = x$levels, as.data.frame(dist),
               as.data.frame(t(std.error)), as.data.frame(do.call("cbind", ci)))
  names(out) <- gsub("theta0", "control.pmf", gsub("theta1", "trt.pmf", names(out)))
  class(out) <- c("tbl_df", "tbl", "data.frame")
  return(out)
}

#' Tidy a(n) mannwhit object
#'
#' @param x A `mannwhit` object.
#' @param ... Unused, included for generic consistency only.
#' @return A tidy [tibble::tibble()] summarizing the Mann-Whitney statistic estimates.
#'
#' @examples
#' rct <- ordinalrct(state_ordinal ~ A + age, A ~ 1, data = c19.ordinal)
#' tidy(mannwhitney(rct))
#'
#' @export
tidy.mannwhit <- function(x, ...) {
  out <- data.frame(mann.whitney = x$estimates$mann.whitney,
                    std.error = x$estimates$std.error[1, 1],
                    conf.low = x$estimates$ci[1],
                    conf.high = x$estimates$ci[2])
  class(out) <- c("tbl_df", "tbl", "data.frame")
  return(out)
}
