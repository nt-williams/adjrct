
compute_rmst <- function(meta, nuis, estimator) {
  switch(estimator,
         tmle = rmst_tmle(meta, nuis),
         aipw = rmst_aipw(meta, nuis),
         ipw  = rmst_ipw(meta, nuis))
}

rmst_tmle <- function(meta, nuis) {

  id   <- access_meta_var(meta, "id")
  trt  <- access_meta_var(meta, "trt")
  lh   <- trt*nuis$H1 + (1 - trt)*nuis$H0
  gr   <- trt*nuis$R1 + (1 - trt)*nuis$R0
  ind  <- outer(meta$m, 1:meta$k, '<=')
  crit <- TRUE
  res  <- list()
  i    <- 1

  for (j in 1:length(meta$tau)) {

    # iterative process
    while (crit && i <= 20) {
      # EIF components
      s1     <- cumprod_by_id(1 - nuis$H1, id)
      s0     <- cumprod_by_id(1 - nuis$H0, id)
      s1_lag <- prodlag_by_id(1 - nuis$H1, id)
      s0_lag <- prodlag_by_id(1 - nuis$H0, id)
      g1     <- cumprod_by_id(1 - nuis$R1, id)
      g0     <- cumprod_by_id(1 - nuis$R0, id)
      z1     <- compute_z(ind, meta$tau[j], id, s1, nuis$A1, g1)
      z0     <- compute_z(ind, meta$tau[j], id, s0, nuis$A0, g0)
      h1     <- compute_z(ind, meta$tau[j], id, s1_lag, nuis$A1, g1)
      h0     <- compute_z(ind, meta$tau[j], id, s0_lag, nuis$A0, g0)
      h      <-  trt * h1 - (1 - trt) * h0
      M      <- compute_m(s1, s0, nuis$A1, nuis$A0, meta$m, meta$tau[j], id)

      # update H
      eps <-
        tilt_params(data.frame(
          lm = meta$data[["lm"]],
          lh = lh,
          A = access_meta_var(meta, "trt"),
          z1 = z1,
          z0 = z0,
          im = meta$im
        ), "epsilon")

      nuis$H1 <- bound01(plogis(qlogis(nuis$H1) + eps[1] * z1))
      nuis$H0 <- bound01(plogis(qlogis(nuis$H0) + eps[2] * z0))

      # update R
      gamma <- tilt_params(data.frame(
        rm = meta$data[["rm"]],
        gr = gr,
        h = h,
        jm = meta$jm
      ), "gamma")

      nuis$R1 <- bound01(plogis(qlogis(nuis$R1) + gamma * h1))
      nuis$R0 <- bound01(plogis(qlogis(nuis$R0) - gamma * h0))

      # update A
      nu <- tilt_params(data.frame(
        A = access_meta_var(meta, "trt"),
        a1 = nuis$A1,
        M = M,
        m = meta$m
      ), "nu")

      nuis$A1 <- bound01(plogis(qlogis(nuis$A1) + nu * M))
      nuis$A0 <- 1 - nuis$A1

      lh <- trt*nuis$H1 + (1 - trt)*nuis$H0
      gr <- trt*nuis$R1 + (1 - trt)*nuis$R0

      # stop criteria update
      i    <-  i + 1
      crit <- any(abs(c(eps, gamma, nu)) > 1e-3/meta$n^(0.6))
    }

    # compute EIF
    s1 <- cumprod_by_id(1 - nuis$H1, id)
    s0 <- cumprod_by_id(1 - nuis$H0, id)
    g1 <- cumprod_by_id(1 - nuis$R1, id)
    g0 <- cumprod_by_id(1 - nuis$R0, id)
    z1 <- compute_z(ind, meta$tau[j], id, s1, nuis$A1, g1)
    z0 <- compute_z(ind, meta$tau[j], id, s0, nuis$A0, g0)

    res[[j]] <- rmst_eif(meta, "tmle", trt, meta$tau[j], z1, z0, s1, s0, lh, id)
  }

  return(res)

}

rmst_aipw <- function(meta, nuis) {
  id  <- access_meta_var(meta, "id")
  trt <- access_meta_var(meta, "trt")
  lh  <- trt*nuis$H1 + (1 - trt)*nuis$H0
  gr  <- trt*nuis$R1 + (1 - trt)*nuis$R0
  ind <- outer(meta$m, 1:meta$k, '<=')
  s1  <- cumprod_by_id(1 - nuis$H1, id)
  s0  <- cumprod_by_id(1 - nuis$H0, id)
  g1  <- cumprod_by_id(1 - nuis$R1, id)
  g0  <- cumprod_by_id(1 - nuis$R0, id)
  z1  <- compute_z(ind, meta$tau, id, s1, nuis$A1, g1)
  z0  <- compute_z(ind, meta$tau, id, s0, nuis$A0, g0)

  rmst_eif(meta, "aipw", trt, z1, z0, s1, s0, lh, id)
}

rmst_ipw <- function(meta, nuis) {
  id  <- access_meta_var(meta, "id")
  trt <- access_meta_var(meta, "trt")
  ind <- outer(meta$m, 1:meta$k, '<=')
  g1  <- cumprod_by_id(1 - nuis$R1, id)
  g0  <- cumprod_by_id(1 - nuis$R0, id)
  z1  <- -rowSums(ind[, 1:(meta$tau - 1)]) / bound(nuis$A1 * unlist(g1))
  z0  <- -rowSums(ind[, 1:(meta$tau - 1)]) / bound(nuis$A0 * unlist(g0))

  rmst_eif(meta, "ipw", trt, z1, z0, NULL, NULL, NULL, id)
}
