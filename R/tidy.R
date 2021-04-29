#' @importFrom generics tidy
#' @export
generics::tidy

#' @export
tidy.rmst <- function(x, ...) {
  out <- x$estimates[-which(names(x$estimates) %in% c("mbcv_theta", "mbcv_treatment", "mbcv_control"))]
  out <- do.call(
    "rbind",
    lapply(out, function(x) as.data.frame(x[-which(names(x) %in% c("eif", "eif1", "eif0"))]))
  )
  out <- cbind(horizon = x$horizon, out)
  class(out) <- c("tbl_df", "tbl", "data.frame")
  return(out)
}

#' @export
tidy.survprob <- tidy.rmst

#' @export
tidy.lor <- function(x, ...) {
  tmp <- x$estimates[-which(names(x$estimates) %in% c("eif", "lo1.ci", "lo0.ci", "ci"))]
  cis <- lapply(x$estimates[c("lo1.ci", "lo0.ci", "ci")], function(x) {
    names(x) <- c("conf.low", "conf.high")
    as.list(x)
  })
  out <- as.data.frame(c(tmp, as.list(unlist(cis))))
  class(out) <- c("tbl_df", "tbl", "data.frame")
  return(out)
}

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
  class(out) <- c("tbl_df", "tbl", "data.frame")
  return(out)
}

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
  class(out) <- c("tbl_df", "tbl", "data.frame")
  return(out)
}

#' @export
tidy.mannwhit <- function(x, ...) {
  out <- data.frame(theta = x$estimates$theta,
                    std.error = x$estimates$std.error[1, 1],
                    conf.low = x$estimates$ci[1],
                    conf.high = x$estimates$ci[2])
  class(out) <- c("tbl_df", "tbl", "data.frame")
  return(out)
}
