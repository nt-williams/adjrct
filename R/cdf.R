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

  crit <- FALSE
  iter <- 1
  ind <- outer(as.numeric(meta$ordinal_data$kl), 1:(K - 1), '<=')
  psin <- Dn <- list()

  while (!crit && iter <= 10) {
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

    H_on <- bound01(plogis(qlogis(H_on) + as.vector(cbind(Z1, 0 * Z0) %*% eps)))
    H_off <- bound01(plogis(qlogis(H_off) + as.vector(cbind(0 * Z1, Z0) %*% eps)))

    Qkn <- trtl*H_on + (1 - trtl)*H_off

    tmpF1 <- matrix(H_on, nrow = m, ncol = K - 1)
    tmpF0 <- matrix(H_off, nrow = m, ncol = K - 1)
    tmpDnY <- Reduce('+', split(as.data.frame(R * Z * (Yl - cbind(tmpF1, tmpF0))), meta$ordinal_data$kl))

    iter <-  iter + 1
    crit <- all(colMeans(tmpDnY) < (vapply(tmpDnY, sd, FUN.VALUE = 1) / (sqrt(n) * log(n)))*0.001)
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

  dist <- compute_theta(H_on, H_off, K, id)
  std.error <- sqrt(sapply(Dnl, function(x) diag(var(x))) / n)

  mbcv1 <- simul::simul(as.list(dist["theta1", ]),
                        apply(sapply(Dnl, function(x) x[, 1]), 2,
                              function(x) x,
                              simplify = FALSE), meta$nobs)

  mbcv0 <- simul::simul(as.list(dist["theta0", ]),
                        apply(sapply(Dnl, function(x) x[, 2]), 2,
                              function(x) x,
                              simplify = FALSE), meta$nobs)

  ci <- list(theta1 = dist_ci(dist[1, ], std.error[1, ]),
             theta0 = dist_ci(dist[2, ], std.error[2, ]),
             unif1 = dist_ci(dist[1, ], std.error[1, ], mbcv1),
             unif0 = dist_ci(dist[2, ], std.error[2, ], mbcv0))

  list(dist = dist,
       std.error = std.error,
       ci = ci,
       eif = list(theta1 = sapply(Dnl, function(x) x[, 1]),
                  theta0 = sapply(Dnl, function(x) x[, 2])),
       mbcv = c("theta1" = mbcv1, "theta0" = mbcv0))
}
