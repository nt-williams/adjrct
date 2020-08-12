
prodlag <- function(x) {
  cumprod(c(1, x[-length(x)]))
}

cumprod_by_id <- function(x, id) {
  tapply(x, id, cumprod, simplify = FALSE)
}

prodlag_by_id <- function(x, id) {
  tapply(x, id, prodlag, simplify = FALSE)
}

sum_by_id <- function(x, id) {
  tapply(x, id, sum)
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

format_digits <- function(x, n) {
  format(round(as.double(x), digits = n), nsmall = n, trim = TRUE)
}

do_rbind <- function(bind, id) {
  do.call("rbind", bind[id])
}
