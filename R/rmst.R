
compute_rmst <- function(meta) {
  switch(meta$estimator,
         tmle = rmst_tmle(meta),
         aipw = rmst_aipw(meta),
         ipw  = rmst_ipw(meta))
}

rmst_tmle <- function(meta) {

  nobs      <- meta$nobs
  id        <- meta$get_var("id")
  trt       <- meta$get_var("trt")
  ind       <- meta$time_indicator()
  evnt      <- meta$surv_data[["evnt"]]
  cens      <- meta$surv_data[["cens"]]
  risk_evnt <- meta$risk_evnt
  risk_cens <- meta$risk_cens
  all_time  <- meta$all_time
  res       <- listenv::listenv()

  for (j in 1:length(meta$horizon)) {
    res[[j]] %<-% {
      time <- meta$horizon[j]
      aux <- Auxiliary$new(meta$nuisance, meta$horizon[j])
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

      rmst_eif("tmle", meta, aux)
    }
  }
  compute_simulband(as.list(res), nobs)
}

rmst_aipw <- function(meta) {

  id  <- meta$get_var("id")
  trt <- meta$get_var("trt")
  ind <- meta$time_indicator()
  res <- listenv::listenv()

  for (j in 1:length(meta$horizon)) {
    res[[j]] %<-% {
      aux  <- Auxiliary$
        new(meta$nuisance, meta$horizon[j])$
        compute_LHGR(trt)$
        compute_S(id)$
        compute_SL(id)$
        compute_G(id)$
        compute_Z_rmst(ind, id)

      rmst_eif("aipw", meta, aux)
    }
  }
  compute_simulband(as.list(res), meta$nobs)
}

# rmst_ipw <- function(meta, nuis) {
#   id  <- access_meta_var(meta, "id")
#   trt <- access_meta_var(meta, "trt")
#   ind <- outer(meta$m, 1:meta$k, '<=')
#   g1  <- cumprod_by_id(1 - nuis$cens_on, id)
#   g0  <- cumprod_by_id(1 - nuis$cens_off, id)
#   z1  <- -rowSums(ind[, 1:(meta$tau - 1)]) / bound(nuis$treat_on * unlist(g1))
#   z0  <- -rowSums(ind[, 1:(meta$tau - 1)]) / bound(nuis$treat_off * unlist(g0))
#
#   rmst_eif(meta, "ipw", trt, z1, z0, NULL, NULL, NULL, id)
# }
