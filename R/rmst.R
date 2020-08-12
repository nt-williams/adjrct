
compute_rmst <- function(metadat) {
  switch(metadat$estimator,
         tmle = rmst_tmle(metadat),
         aipw = rmst_aipw(metadat),
         ipw  = rmst_ipw(metadat))
}

rmst_tmle <- function(metadat) {

  id   <- access_meta_var(metadat, "id")
  trt  <- access_meta_var(metadat, "trt")
  ind  <- outer(metadat$m, 1:metadat$k, '<=')
  res  <- list()

  for (j in 1:length(metadat$horizon)) {
    updated_nuis <- metadat$nuisance
    time         <- metadat$horizon[j]
    lh           <- trt*updated_nuis$hazard_on + (1 - trt)*updated_nuis$hazard_off
    gr           <- trt*updated_nuis$cens_on + (1 - trt)*updated_nuis$cens_off
    crit         <- TRUE
    i            <- 1
    while (crit && i <= 20) {
      s1     <- cumprod_by_id(1 - updated_nuis$hazard_on, id)
      s0     <- cumprod_by_id(1 - updated_nuis$hazard_off, id)
      s1_lag <- prodlag_by_id(1 - updated_nuis$hazard_on, id)
      s0_lag <- prodlag_by_id(1 - updated_nuis$hazard_off, id)
      g1     <- cumprod_by_id(1 - updated_nuis$cens_on, id)
      g0     <- cumprod_by_id(1 - updated_nuis$cens_off, id)
      z1     <- compute_z(ind, time, id, s1, updated_nuis$treat_on, g1)
      z0     <- compute_z(ind, time, id, s0, updated_nuis$treat_off, g0)
      h1     <- compute_z(ind, time, id, s1_lag, updated_nuis$treat_on, g1)
      h0     <- compute_z(ind, time, id, s0_lag, updated_nuis$treat_off, g0)
      h      <- trt * h1 - (1 - trt) * h0
      M      <- compute_m(s1, s0, updated_nuis$treat_on, updated_nuis$treat_off, metadat$m, time, id)

      eps <-
        tilt_params(data.frame(lm = metadat$data[["lm"]],
                               lh = lh,
                               A = trt,
                               z1 = z1,
                               z0 = z0,
                               im = metadat$im),
                    "epsilon")

      updated_nuis$hazard_on  <- bound01(plogis(qlogis(updated_nuis$hazard_on) + eps[1] * z1))
      updated_nuis$hazard_off <- bound01(plogis(qlogis(updated_nuis$hazard_off) + eps[2] * z0))

      gamma <- tilt_params(data.frame(rm = metadat$data[["rm"]],
                                      gr = gr,
                                      h = h,
                                      jm = metadat$jm),
                           "gamma")

      updated_nuis$cens_on  <- bound01(plogis(qlogis(updated_nuis$cens_on) + gamma * h1))
      updated_nuis$cens_off <- bound01(plogis(qlogis(updated_nuis$cens_off) - gamma * h0))

      nu <- tilt_params(data.frame(A = trt,
                                   a1 = updated_nuis$treat_on,
                                   M = M,
                                   m = metadat$m),
                        "nu")

      updated_nuis$treat_on  <- bound01(plogis(qlogis(updated_nuis$treat_on) + nu * M))
      updated_nuis$treat_off <- 1 - updated_nuis$treat_on

      lh   <- trt*updated_nuis$hazard_on + (1 - trt)*updated_nuis$hazard_off
      gr   <- trt*updated_nuis$cens_on + (1 - trt)*updated_nuis$cens_off
      i    <-  i + 1
      crit <- any(abs(c(eps, gamma, nu)) > 1e-3/metadat$n^(0.6))
    }

    s1 <- cumprod_by_id(1 - updated_nuis$hazard_on, id)
    s0 <- cumprod_by_id(1 - updated_nuis$hazard_off, id)
    g1 <- cumprod_by_id(1 - updated_nuis$cens_on, id)
    g0 <- cumprod_by_id(1 - updated_nuis$cens_off, id)
    z1 <- compute_z(ind, time, id, s1, updated_nuis$treat_on, g1)
    z0 <- compute_z(ind, time, id, s0, updated_nuis$treat_off, g0)

    res[[j]] <- rmst_eif(metadat, "tmle", trt, time, z1, z0, s1, s0, lh, id)
  }
  res <- compute_simulband(res, metadat$n)
  return(res)
}

rmst_aipw <- function(meta, nuis) {
  id  <- access_meta_var(meta, "id")
  trt <- access_meta_var(meta, "trt")
  lh  <- trt*nuis$hazard_on + (1 - trt)*nuis$hazard_off
  gr  <- trt*nuis$cens_on + (1 - trt)*nuis$cens_off
  ind <- outer(meta$m, 1:meta$k, '<=')
  s1  <- cumprod_by_id(1 - nuis$hazard_on, id)
  s0  <- cumprod_by_id(1 - nuis$hazard_off, id)
  g1  <- cumprod_by_id(1 - nuis$cens_on, id)
  g0  <- cumprod_by_id(1 - nuis$cens_off, id)
  z1  <- compute_z(ind, meta$tau, id, s1, nuis$treat_on, g1)
  z0  <- compute_z(ind, meta$tau, id, s0, nuis$treat_off, g0)

  rmst_eif(meta, "aipw", trt, z1, z0, s1, s0, lh, id)
}

rmst_ipw <- function(meta, nuis) {
  id  <- access_meta_var(meta, "id")
  trt <- access_meta_var(meta, "trt")
  ind <- outer(meta$m, 1:meta$k, '<=')
  g1  <- cumprod_by_id(1 - nuis$cens_on, id)
  g0  <- cumprod_by_id(1 - nuis$cens_off, id)
  z1  <- -rowSums(ind[, 1:(meta$tau - 1)]) / bound(nuis$treat_on * unlist(g1))
  z0  <- -rowSums(ind[, 1:(meta$tau - 1)]) / bound(nuis$treat_off * unlist(g0))

  rmst_eif(meta, "ipw", trt, z1, z0, NULL, NULL, NULL, id)
}
