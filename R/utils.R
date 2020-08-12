
prodlag <- function(x) {
  cumprod(c(1, x[-length(x)]))
}

cumprod_by_id <- function(x, id) {
  tapply(x, id, cumprod, simplify = FALSE)
}

prodlag_by_id <- function(x, id) {
  tapply(x, id, prodlag, simplify = FALSE)
}

compute_z <- function(ind, time, id, S, A, R) {
  -rowSums(matrix((ind * do.call("rbind", S[id]))[, 1:(time - 1)],
                  nrow = length(id),
                  ncol = time - 1)) /
    bound(unlist(S) * A[id] * unlist(R))
}

compute_z_survprob <- function(ind, id, time, S, A, R) {
  -(ind * do.call("rbind", S[id]))[, time] /
    bound(unlist(S) * A[id] * unlist(R))
}

compute_h_survprob <- function(id, time, S, A, R) {
  -do.call("rbind", S[id])[, time] /
    bound(unlist(S) * A[id] * unlist(R))
}

sum_by_id <- function(x, id) {
  tapply(x, id, sum)
}

compute_m <- function(S1, S0, A1, A0, m, time, id) {
  (sum_by_id(unlist(S1) / bound(A1[id]) * (m <= time - 1), id) +
    sum_by_id(unlist(S0) / bound(A0[id]) * (m <= time - 1), id))[id]
}

compute_m_survprob <- function(S1, S0, A1, A0, m, time, id) {
  (sum_by_id(unlist(S1) / bound(A1[id]) * (m <= time), id) +
     sum_by_id(unlist(S0) / bound(A0[id]) * (m <= time), id))[id]
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
