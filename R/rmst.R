
compute_rmst <- function(metadat) {
  switch(metadat$estimator,
         tmle = rmst_tmle(metadat),
         aipw = rmst_aipw(metadat),
         ipw  = rmst_ipw(metadat))
}

rmst_tmle <- function(metadat) {

  nobs      <- metadat$nobs
  id        <- metadat$get_var("id")
  trt       <- metadat$get_var("trt")
  ind       <- metadat$time_indicator()
  evnt      <- metadat$surv_data[["evnt"]]
  cens      <- metadat$surv_data[["cens"]]
  risk_evnt <- metadat$risk_evnt
  risk_cens <- metadat$risk_cens
  all_time  <- metadat$all_time
  res       <- listenv::listenv()

  for (j in 1:length(metadat$horizon)) {
    res[[j]] %<-% {
      time <- metadat$horizon[j]
      aux <- Auxiliary$new(metadat$nuisance, metadat$horizon[j])
      while (aux$crit && aux$iter <= 20) {
        aux$
          compute_LHGR(trt)$
          compute_S(id)$
          compute_SL(id)$
          compute_G(id)$
          compute_Z_rmst(ind, id)$
          compute_H_rmst(ind, trt, id)$
          compute_M_rmst(id, all_time)$
          tilt_eps(trt, evnt, risk_evnt)$
          tilt_gamma(cens, risk_cens)$
          tilt_nu(trt, all_time)$
          update_crit(nobs)
      }
      aux$
        compute_S(id)$
        compute_G(id)$
        compute_Z_rmst(ind, id)$
        compute_LHGR(trt)

      rmst_eif("tmle", metadat, aux)
    }
  }

  res <- compute_simulband(as.list(res), nobs)
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
