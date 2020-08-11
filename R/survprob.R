
compute_survprob <- function(meta, nuis, estimator) {
  switch(estimator,
         tmle = survprob_tmle(meta, nuis))
}

survprob_tmle <- function(meta, nuis) {

  id   <- access_meta_var(meta, "id")
  trt  <- access_meta_var(meta, "trt")
  ind  <- outer(meta$m, 1:meta$k, '<=')
  res  <- list()

  for (j in 1:length(meta$tau)) {

    time         <- meta$tau[j]
    updated_nuis <- nuis
    lh           <- trt*updated_nuis$hazard_on + (1 - trt)*updated_nuis$hazard_off
    gr           <- trt*updated_nuis$cens_on + (1 - trt)*updated_nuis$cens_off
    crit         <- TRUE
    i            <- 1

    # iterative process
    while (crit && i <= 20) {
      # EIF components
      s1     <- cumprod_by_id(1 - updated_nuis$hazard_on, id)
      s0     <- cumprod_by_id(1 - updated_nuis$hazard_off, id)
      s1_lag <- prodlag_by_id(1 - updated_nuis$hazard_on, id)
      s0_lag <- prodlag_by_id(1 - updated_nuis$hazard_off, id)
      g1     <- cumprod_by_id(1 - updated_nuis$cens_on, id)
      g0     <- cumprod_by_id(1 - updated_nuis$cens_off, id)
      z1     <- compute_z_survprob(ind, id, time, s1, updated_nuis$treat_on, g1)
      z0     <- compute_z_survprob(ind, id, time, s0, updated_nuis$treat_off, g0)
      h1     <- compute_h_survprob(id, time, s1_lag, updated_nuis$treat_on, g1)
      h0     <- compute_h_survprob(id, time, s0_lag, updated_nuis$treat_off, g0)
      h      <-  trt * h1 - (1 - trt) * h0
      M      <- compute_m_survprob(s1, s0, updated_nuis$treat_on, updated_nuis$treat_off, meta$m, time, id)

      # update H
      eps <-
        tilt_params(data.frame(lm = meta$data[["lm"]],
                               lh = lh,
                               A = access_meta_var(meta, "trt"),
                               z1 = z1,
                               z0 = z0,
                               im = meta$im),
                    "epsilon")

      updated_nuis$hazard_on  <- bound01(plogis(qlogis(updated_nuis$hazard_on) + eps[1] * z1))
      updated_nuis$hazard_off <- bound01(plogis(qlogis(updated_nuis$hazard_off) + eps[2] * z0))

      # update R
      gamma <- tilt_params(data.frame(rm = meta$data[["rm"]],
                                      gr = gr,
                                      h = h,
                                      jm = meta$jm),
                           "gamma")

      updated_nuis$cens_on  <- bound01(plogis(qlogis(updated_nuis$cens_on) + gamma * h1))
      updated_nuis$cens_off <- bound01(plogis(qlogis(updated_nuis$cens_off) - gamma * h0))

      # update A
      nu <- tilt_params(data.frame(A = access_meta_var(meta, "trt"),
                                   a1 = updated_nuis$treat_on,
                                   M = M,
                                   m = meta$m),
                        "nu")

      updated_nuis$treat_on  <- bound01(plogis(qlogis(updated_nuis$treat_on) + nu * M))
      updated_nuis$treat_off <- 1 - updated_nuis$treat_on

      lh <- trt*updated_nuis$hazard_on + (1 - trt)*updated_nuis$hazard_off
      gr <- trt*updated_nuis$cens_on + (1 - trt)*updated_nuis$cens_off

      # stop criteria update
      i    <-  i + 1
      crit <- any(abs(c(eps, gamma, nu)) > 1e-3/meta$n^(0.6))
    }

    # compute EIF
    s1 <- cumprod_by_id(1 - updated_nuis$hazard_on, id)
    s0 <- cumprod_by_id(1 - updated_nuis$hazard_off, id)
    g1 <- cumprod_by_id(1 - updated_nuis$cens_on, id)
    g0 <- cumprod_by_id(1 - updated_nuis$cens_off, id)
    z1 <- compute_z_survprob(ind, id, time, s1, updated_nuis$treat_on, g1)
    z0 <- compute_z_survprob(ind, id, time, s0, updated_nuis$treat_off, g0)

    res[[j]] <- survprob_eif(meta, "tmle", trt, time, z1, z0, s1, s0, lh, id)

  }
  # multiplier bootstrap
  res <- compute_simulband(res, meta$n)

  return(res)
}
