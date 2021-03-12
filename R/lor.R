compute_theta <- function(H_on, H_off, K, id) {
  tmp1 <- tapply(1 - H_on, id, cumprod, simplify = FALSE)
  tmp0 <- tapply(1 - H_off, id, cumprod, simplify = FALSE)

  if (K > 2) {
    theta1 <- 1 - rowMeans(sapply(tmp1, I))
    theta0 <- 1 - rowMeans(sapply(tmp0, I))
  }
  if (K == 2) {
    theta1 <- 1 - mean(sapply(tmp1, I))
    theta0 <- 1 - mean(sapply(tmp0, I))
  }
  rbind(theta1 = theta1, theta0 = theta0)
}

compute_lor <- function(x) {
  switch(x$estimator,
         tmle = lor_tmle(x),
         unadjusted = lor_unadjusted(x))
}

lor_tmle <- function(meta) {
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

  Qkn <- trtl*H_on + (1 - trtl)*H_off
  gAn <- trtl*trt_on + (1 - trtl)*trt_off

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

    theta <- compute_theta(H_on, H_off, K, id)

    Z1 <- ind * prodk1 / (prodj1 * trt_on[id])
    Z0 <- ind * prodk0 / (prodj0 * trt_off[id])
    Z <- cbind(trtl * Z1, (1 - trtl) * Z0)

    M1 <- do.call(rbind, tapply(1 - H_on, id, cumprod, simplify = FALSE))
    M0 <- do.call(rbind, tapply(1 - H_off, id, cumprod, simplify = FALSE))

    H1 <- (M1 - theta[1,]) / (trt_on)
    H0 <- (M0 - theta[2,]) / (trt_off)
    H1 <- colSums(t(H1) / (theta[1, ] * (1 - theta[1, ])))
    H0 <- colSums(t(H0) / (theta[2, ] * (1 - theta[2,])))
    H  <- trt * H1 - (1 - trt) * H0

    M1 <- (M1 - theta[1,]) / trt_on
    M0 <- (M0 - theta[2,]) / trt_off
    M1 <- colSums(t(M1) / (theta[1,] * (1 - theta[1, ])))
    M0 <- colSums(t(M0) / (theta[2,] * (1 - theta[2, ])))
    M <- M1 + M0

    eps <- coef(glm(Yl ~ 0 + offset(qlogis(Qkn)) + Z, family = binomial(), subset = R == 1, start = rep(0, ncol(Z))))
    nu <- coef(glm(trt ~ 0 + offset(qlogis(trt_on)) + M, family = binomial(), start = 0))

    eps[is.na(eps)] <- 0
    nu[is.na(nu)] <- 0

    H_on <- bound01(plogis(qlogis(H_on) + as.vector(cbind(Z1, 0 * Z0) %*% eps)))
    H_off <- bound01(plogis(qlogis(H_off) + as.vector(cbind(0 * Z1, Z0) %*% eps)))

    trt_on <- bound01(plogis(qlogis(trt_on) + M * nu))
    trt_off <- 1 - trt_on

    Qkn <- trtl*H_on + (1 - trtl)*H_off
    gAn <- trt*trt_on + (1 - trt)*trt_off

    tmpZ1 <- colSums(t(Z1) / (theta[1,] * (1 - theta[1,])))
    tmpZ0 <- colSums(t(Z0) / (theta[2,] * (1 - theta[2,])))
    tmpZ <- trtl * tmpZ1 - (1 - trtl) * tmpZ0
    tmpDnY <- tapply(R * tmpZ * (Yl - Qkn), id, mean)

    iter <- iter + 1
    crit <- mean(tmpDnY) < (sd(tmpDnY) / (sqrt(n) * log(n)))*0.001
  }

  Z1 <- colSums(t(Z1) / (theta[1,] * (1 - theta[1,])))
  Z0 <- colSums(t(Z0) / (theta[2,] * (1 - theta[2,])))
  Z <- trtl * Z1 - (1 - trtl) * Z0
  DnY <- tapply(R * Z * (Yl - Qkn), id, mean)
  DnY1 <- tapply(R * trtl*Z1 * (Yl - Qkn), id, mean)
  DnY0 <- tapply(R * (1 - trtl)*Z0 * (Yl - Qkn), id, mean)

  tmp1 <- colMeans(do.call(cbind, tapply(1 - H_on, id, cumprod, simplify = FALSE))
                   / (theta[1, ] * (1 - theta[1, ])))
  tmp0 <- colMeans(do.call(cbind, tapply(1 - H_off, id, cumprod, simplify = FALSE))
                   / (theta[2, ] * (1 - theta[2, ])))

  eif <- DnY - tmp1 + tmp0
  eif1 <- DnY1 - tmp1
  eif0 <- DnY0 - tmp0
  cdf <- compute_theta(H_on, H_off, K, id)

  lor <- LOR(cdf)
  std.error <- sqrt(var(eif) / n)
  lo1.std.error <- sqrt(var(eif1) / n)
  lo0.std.error <- sqrt(var(eif0) / n)

  list(lo1 = lor$arm1,
       lo1.std.error = lo1.std.error,
       lo1.ci = c(lor$arm1 - qnorm(0.975)*lo1.std.error, lor$arm1 + qnorm(0.975)*lo1.std.error),
       lo0 = lor$arm0,
       lo0.std.error = lo0.std.error,
       lo0.ci = c(lor$arm0 - qnorm(0.975)*lo0.std.error, lor$arm0 + qnorm(0.975)*lo0.std.error),
       lor = lor$theta,
       std.error = std.error,
       ci = c(lor$theta - qnorm(0.975)*std.error, lor$theta + qnorm(0.975)*std.error),
       eif = eif)
}

tcs <- function(x) {
  out <- vector("numeric", length(x))
  for (i in 1:length(x)) {
    out[i] <- sum(x[1:i], na.rm = T)
  }
  out
}

LOR <- function(cdf) {
  theta1 <- mean(log(cdf[1, ] / (1 - cdf[1, ])))
  theta0 <- mean(log((cdf[2, ]) / (1 - cdf[2, ])))
  list(arm1 = theta1, arm0 = theta0, theta = theta1 - theta0)
}

lor_unadjusted <- function(meta) {
  tab <- table(meta$data[[meta$Y]], meta$data[[meta$trt]])
  tabcumsum <- apply(tab, 2, tcs)
  cdf <- cbind(tabcumsum[1:(nrow(tabcumsum) - 1), 2] /
                 sum(tab[, 2], na.rm = T), tabcumsum[1:(nrow(tabcumsum) - 1), 1] /
                 sum(tab[, 1], na.rm = T))
  list(lor = LOR(t(cdf)),
       std.error = NULL,
       ci = NULL)
}
