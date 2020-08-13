
compute_survprob <- function(meta, nuis, estimator) {
  switch(meta$estimator,
         tmle = survprob_tmle(meta, nuis),
         aipw = survprob_aipw(meta, nuis))
}

survprob_tmle <- function(meta, nuis) {

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
          compute_Z_survprob(ind, id)$
          compute_H_survprob(id, trt)$
          compute_M_survprob(id, all_time)$
          tilt_eps(trt, evnt, risk_evnt)$
          tilt_gamma(cens, risk_cens)$
          tilt_nu(trt, all_time)$
          update_crit(nobs)
      }
      aux$
        compute_S(id)$
        compute_G(id)$
        compute_Z_survprob(ind, id)$
        compute_LHGR(trt)

      survprob_eif("tmle", meta, aux)
    }
  }
  compute_simulband(as.list(res), nobs)
}

survprob_aipw <- function(meta, nuis) {

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
        compute_Z_survprob(ind, id)

      survprob_eif("aipw", meta, aux)
    }
  }
  compute_simulband(as.list(res), meta$nobs)
}
