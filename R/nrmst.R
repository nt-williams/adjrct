compute_rmst <- function(x) {
  switch(x$estimator,
         tmle = rmst_tmle(x))
}

rmst_tmle <- function(meta) {
  nobs <- meta$nobs
  id <- meta$surv_data[["survrctId"]]
  trt <- meta$get_var("trt")
  ind <- meta$time_indicator()
  evnt <- meta$surv_data[["evnt"]]
  cens <- meta$surv_data[["cens"]]
  risk_evnt <- meta$risk_evnt
  risk_cens <- meta$risk_cens
  all_time <- meta$all_time
  nuis <- meta$nuisance
  res <- listenv::listenv()

  gA1 <- nuis$trt_on
  gA0 <- 1 - gA1

  h1 <- nuis$hzrd_on
  h0 <- nuis$hzrd_off
  gR1 <- nuis$cens_on
  gR0 <- nuis$cens_off

  h <- trt*h1 + (1 - trt)*h0
  gR <- trt*gR1 + (1 - trt)*gR0

  for (j in 1:length(meta$horizon)) {
    res[[j]] %<-% {
      tau <- meta$horizon[j]
      crit <- TRUE
      iter <- 1
      while (crit && iter <= 10) {
        S1 <- tapply(1 - h1, id, cumprod, simplify = FALSE)
        S0 <- tapply(1 - h0, id, cumprod, simplify = FALSE)
        S1lag <- tapply(1 - h1, id, prodlag, simplify = FALSE)
        S0lag <- tapply(1 - h0, id, prodlag, simplify = FALSE)

        G1 <- tapply(1 - gR1, id, cumprod, simplify = FALSE)
        G0 <- tapply(1 - gR0, id, cumprod, simplify = FALSE)

        St1 <- do.call('rbind', S1[id])
        St0 <- do.call('rbind', S0[id])

        Sm1 <- unlist(S1)
        Sm0 <- unlist(S0)
        Sm1lag <- unlist(S1lag)
        Sm0lag <- unlist(S0lag)

        Gm1 <- unlist(G1)
        Gm0 <- unlist(G0)

        Z1 <- -rowSums((ind * St1)[, 1:(tau - 1), drop = FALSE]) / bound(Sm1 * gA1[id] * Gm1)
        Z0 <- -rowSums((ind * St0)[, 1:(tau - 1), drop = FALSE]) / bound(Sm0 * gA0[id] * Gm0)

        M <- tapply(Sm1 / bound(gA1[id]) * (all_time <= tau - 1), id, sum) +
          tapply(Sm0 / bound(gA0[id]) * (all_time <= tau - 1), id, sum)

        H1 <- -rowSums((ind * St1)[, 1:(tau - 1), drop = FALSE]) / bound(Sm1lag * gA1[id] * Gm1)
        H0 <- -rowSums((ind * St0)[, 1:(tau - 1), drop = FALSE]) / bound(Sm0lag * gA0[id] * Gm0)
        H <- trt * H1 - (1 - trt) * H0

        eps <- coef(glm2::glm2(evnt[risk_evnt == 1] ~ 0 + offset(qlogis(h[risk_evnt == 1])) + I(trt[risk_evnt == 1] * Z1[risk_evnt == 1]) + I((1-trt[risk_evnt == 1]) * Z0[risk_evnt == 1]), family = binomial()))
        gamma <- coef(glm2::glm2(cens[risk_cens == 1] ~ 0 + offset(qlogis(gR[risk_cens == 1])) + H[risk_cens == 1], family = binomial()))
        nu <- coef(glm2::glm2(trt[all_time == 1] ~ 0 + offset(qlogis(gA1)) + M, family = binomial()))

        eps[is.na(eps)] <- 0
        gamma[is.na(gamma)] <- 0
        nu[is.na(nu)] <- 0

        h1 <- bound01(plogis(qlogis(h1) + eps[1] * Z1))
        h0 <- bound01(plogis(qlogis(h0) + eps[2] * Z0))
        gR1 <- bound01(plogis(qlogis(gR1) + gamma * H1))
        gR0 <- bound01(plogis(qlogis(gR0) - gamma * H0))
        gA1 <- bound01(plogis(qlogis(gA1) + nu * M))
        gA0 <- 1 - gA1
        h <- trt*h1 + (1 - trt)*h0
        gR <- trt*gR1 + (1 - trt)*gR0

        iter <- iter + 1
        crit <- any(abs(c(eps, gamma, nu)) > 1e-3/nobs^(0.6))
      }

      S1 <- tapply(1 - h1, id, cumprod, simplify = FALSE)
      S0 <- tapply(1 - h0, id, cumprod, simplify = FALSE)
      G1 <- tapply(1 - gR1, id, cumprod, simplify = FALSE)
      G0 <- tapply(1 - gR0, id, cumprod, simplify = FALSE)
      St1 <- do.call('rbind', S1[id])
      St0 <- do.call('rbind', S0[id])
      Sm1 <- unlist(S1)
      Sm0 <- unlist(S0)
      Gm1 <- unlist(G1)
      Gm0 <- unlist(G0)
      Z1 <- -rowSums((ind * St1)[, 1:(tau - 1), drop = FALSE]) / bound(Sm1 * gA1[id] * Gm1)
      Z0 <- -rowSums((ind * St0)[, 1:(tau - 1), drop = FALSE]) / bound(Sm0 * gA0[id] * Gm0)

      DT1 <- tapply(risk_evnt*trt*Z1*(evnt - h), id, sum)
      DT0 <- tapply(risk_evnt*(1 - trt)*Z0*(evnt - h), id, sum)
      DW1 <- rowSums(St1[all_time == 1, 1:(tau - 1), drop = FALSE])
      DW0 <- rowSums(St0[all_time == 1, 1:(tau - 1), drop = FALSE])

      theta1 <- 1 + mean(DW1)
      theta0 <- 1 + mean(DW0)

      eif1 <- as.vector(DT1 + DW1)
      eif0 <- as.vector(DT0 + DW0)
      eif <- eif1 - eif0

      se1 <- sqrt(var(eif1) / nobs)
      se0 <- sqrt(var(eif0) / nobs)
      se <- sqrt(var(eif) / nobs)

      list(arm1 = theta1,
           eif1 = eif1,
           arm1.std.error = se1,
           arm1.conf.low = theta1 - qnorm(0.975)*se1,
           arm1.conf.high = theta1 + qnorm(0.975)*se1,
           arm0 = theta0,
           eif0 = eif0,
           arm0.std.error = se0,
           arm0.conf.low = theta0 - qnorm(0.975)*se0,
           arm0.conf.high = theta0 + qnorm(0.975)*se0,
           theta = theta1 - theta0,
           eif = eif,
           std.error = se,
           theta.conf.low  = (theta1 - theta0) - qnorm(0.975)*se,
           theta.conf.high = (theta1 - theta0) + qnorm(0.975)*se)
    }
  }
  simul_ci(as.list(res), nobs)
}
