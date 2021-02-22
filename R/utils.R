prodlag <- function(x) {
  cumprod(c(1, x[-length(x)]))
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

format_digits <- function(x, n) {
  format(round(as.double(x), digits = n), nsmall = n, trim = TRUE)
}

get_status <- function(formula) {
  all.vars(formula[[2]])[[2]]
}

get_time <- function(formula) {
  all.vars(formula[[2]])[[1]]
}

get_covar <- function(formula, target) {
  setdiff(all.vars(formula[[3]]), target)
}
