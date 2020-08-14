
compute_rmst <- function(meta) {
  switch(meta$estimator,
         tmle = rmst_tmle(meta),
         aipw = rmst_ee(meta),
         km   = rmst_ee(meta))
}

rmst_tmle <- function(meta) {

  nobs      <- meta$nobs
  id        <- meta$surv_data[["survrctId"]]
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

      rmst_eif(meta, aux)
    }
  }
  compute_simulband(as.list(res), nobs)
}

rmst_ee <- function(meta) {

  id  <- meta$surv_data[["survrctId"]]
  trt <- meta$get_var("trt")
  ind <- meta$time_indicator()
  res <- listenv::listenv()

  for (j in 1:length(meta$horizon)) {
    res[[j]] %<-% {
      aux <- Auxiliary$
        new(meta$nuisance, meta$horizon[j])$
        compute_LHGR(trt)$
        compute_S(id)$
        compute_G(id)$
        compute_Z_rmst(ind, id)

      rmst_eif(meta, aux)
    }
  }
  compute_simulband(as.list(res), meta$nobs)
}
