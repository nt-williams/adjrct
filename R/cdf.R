compute_cdf <- function(x) {
  switch(x$estimator,
         tmle = cdf_tmle(x))
}

cdf_tmle <- function(meta) {
  trt <- meta$data[[meta$trt]]
  trtl <- meta$ordinal_data[["A"]]
  Y <- meta$data[[meta$Y]]
  K <- meta$K
  n <- meta$nobs
  id <- meta$id
  R <- meta$R
  Yl <- meta$ordinal_data$Y
  nuis <- meta$nuisance
  H_on <- nuis$H_on
  H_off <- nuis$H_off
  trt_on <- nuis$trt_on
  trt_off <- nuis$trt_off
  m <- length(id)

  Qkn <- trtl*H_on + (1 - trtl)*H_off

  crit <- TRUE
  iter <- 1
  ind <- outer(as.numeric(meta$ordinal_data$kl), 1:(K - 1), '<=')
  psin <- Dn <- list()

  while (crit && iter <= 100) {
    tmp1 <- tapply(1 - H_on, id, cumprod, simplify = FALSE)
    tmp0 <- tapply(1 - H_off, id, cumprod, simplify = FALSE)
    prodk1 <- do.call('rbind', tmp1[id])
    prodk0 <- do.call('rbind', tmp0[id])
    prodj1 <- unlist(tmp1)
    prodj0 <- unlist(tmp0)

    Z1 <- ind * prodk1 / (prodj1 * trt_on[id])
    Z0 <- ind * prodk0 / (prodj0 * trt_off[id])
    Z <- cbind(trtl * Z1, (1 - trtl) * Z0)

    eps <- coef(glm(Yl ~ 0 + offset(qlogis(Qkn)) + Z, family = binomial(), subset = R == 1))

    H_on <- bound(plogis(qlogis(H_on) + as.vector(cbind(Z1, 0 * Z0) %*% eps)))
    H_off <- bound(plogis(qlogis(H_off) + as.vector(cbind(0 * Z1, Z0) %*% eps)))

    Qkn <- trtl*H_on + (1 - trtl)*H_off

    iter <-  iter + 1
    crit <- any(abs(eps) > 1e-3/n^(0.6))
  }

  tmp1 <- matrix(H_on, nrow = m, ncol = K - 1)
  tmp0 <- matrix(H_off, nrow = m, ncol = K - 1)
  DnY <- Reduce('+', split(as.data.frame(R * Z * (Yl - cbind(tmp1, tmp0))), meta$ordinal_data$kl))
  tmp1 <- do.call(rbind, tapply(1 - H_on, id, cumprod, simplify = FALSE))
  tmp0 <- do.call(rbind, tapply(1 - H_off, id, cumprod, simplify = FALSE))

  Dnl <- list()
  for (k in 1:(K - 1)) {
    Dnl[[k]] <- matrix(0, nrow = n, ncol = 2)
    Dnl[[k]][,1] <- DnY[, k] - tmp1[, k]
    Dnl[[k]][,2] <- DnY[, k + K - 1] - tmp0[, k]
  }

  std.error <- sqrt(sapply(Dnl, function(x) diag(var(x))) / n)

  list(dist = compute_theta(H_on, H_off, K, id),
       std.error = std.error,
       eif = list(theta1 = sapply(Dnl, function(x) x[, 1]),
                  theta0 = sapply(Dnl, function(x) x[, 2])))
}
