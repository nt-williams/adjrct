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
tidy.survprob <- function(x, ...) {
  tidy.rmst(x)
}

tidy.lor <- function(x, ...) {

}

tidy.cdf <- function(x, ...) {

}

tidy.pmf <- function(x, ...) {

}

tidy.mannwhit <- function(x, ...) {

}
