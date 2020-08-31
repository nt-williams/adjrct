
Auxiliary <- R6::R6Class(
  "Auxiliary",
  public = list(
    time = NULL,
    crit = TRUE,
    iter = 1,
    nuis = list(),
    LH = NULL,
    GR = NULL,
    S1 = NULL,
    S0 = NULL,
    SL1 = NULL,
    SL0 = NULL,
    G1 = NULL,
    G0 = NULL,
    Z1 = NULL,
    Z0 = NULL,
    H1 = NULL,
    H0 = NULL,
    H = NULL,
    M = NULL,
    eps = NULL,
    gamma = NULL,
    nu = NULL,
    initialize = function(nuis, time) {
      self$time <- time
      self$nuis <- nuis
    },
    compute_LHGR = function(trt) {
      self$LH <- trt*self$nuis$hzrd_on + (1 - trt)*self$nuis$hzrd_off
      self$GR <- trt*self$nuis$cens_on + (1 - trt)*self$nuis$cens_off
      invisible(self)
    },
    compute_S = function(id) {
      self$S1 <- cumprod_by_id(1 - self$nuis$hzrd_on, id)
      self$S0 <- cumprod_by_id(1 - self$nuis$hzrd_off, id)
      invisible(self)
    },
    compute_SL = function(id) {
      self$SL1 <- prodlag_by_id(1 - self$nuis$hzrd_on, id)
      self$SL0 <- prodlag_by_id(1 - self$nuis$hzrd_off, id)
      invisible(self)
    },
    compute_G = function(id) {
      self$G1 <- cumprod_by_id(1 - self$nuis$cens_on, id)
      self$G0 <- cumprod_by_id(1 - self$nuis$cens_off, id)
      invisible(self)
    },
    compute_Z_rmst = function(time_ind, id) {
      self$Z1 <- -rowSums(matrix((time_ind * do_rbind(self$S1, id))[, 1:(self$time - 1)],
                                 nrow = length(id),
                                 ncol = self$time - 1)) /
        bound(unlist(self$S1) * self$nuis$trt_on[id] * unlist(self$G1))
      self$Z0 <- -rowSums(matrix((time_ind * do_rbind(self$S0, id))[, 1:(self$time - 1)],
                                 nrow = length(id),
                                 ncol = self$time - 1)) /
        bound(unlist(self$S0) * self$nuis$trt_off[id] * unlist(self$G0))
      invisible(self)
    },
    compute_H_rmst = function(time_ind, trt, id) {
      self$H1 <- -rowSums(matrix((time_ind * do_rbind(self$SL1, id))[, 1:(self$time - 1)],
                                 nrow = length(id),
                                 ncol = self$time - 1)) /
        bound(unlist(self$SL1) * self$nuis$trt_on[id] * unlist(self$G1))
      self$H0 <- -rowSums(matrix((time_ind * do_rbind(self$SL0, id))[, 1:(self$time - 1)],
                                 nrow = length(id),
                                 ncol = self$time - 1)) /
        bound(unlist(self$SL0) * self$nuis$trt_off[id] * unlist(self$G0))
      self$H <- trt * self$H1 - (1 - trt) * self$H0
      invisible(self)
    },
    compute_Z_survprob = function(time_ind, id) {
      self$Z1 <- -(time_ind * do_rbind(self$S1, id))[, self$time] /
        bound(unlist(self$S1) * self$nuis$trt_on[id] * unlist(self$G1))
      self$Z0 <- -(time_ind * do_rbind(self$S0, id))[, self$time] /
        bound(unlist(self$S0) * self$nuis$trt_off[id] * unlist(self$G0))
      invisible(self)
    },
    compute_H_survprob = function(id, trt) {
      self$H1 <- -do_rbind(self$S1, id)[, self$time] /
        bound(unlist(self$S1) * self$nuis$trt_on[id] * unlist(self$G1))
      self$H0 <- -do_rbind(self$S0, id)[, self$time] /
        bound(unlist(self$S0) * self$nuis$trt_off[id] * unlist(self$G0))
      self$H <- trt * self$H1 - (1 - trt) * self$H0
      invisible(self)
    },
    compute_M_rmst = function(id, current_time) {
      self$M <- (sum_by_id(unlist(self$S1) / bound(self$nuis$trt_on[id]) * (current_time <= self$time - 1), id) +
                   sum_by_id(unlist(self$S0) / bound(self$nuis$trt_off[id]) * (current_time <= self$time - 1), id))[id]
      invisible(self)
    },
    compute_M_survprob = function(id, current_time) {
      self$M <- (sum_by_id(unlist(self$S1) / bound(self$nuis$trt_on[id]) * (current_time <= self$time), id) +
                   sum_by_id(unlist(self$S0) / bound(self$nuis$trt_off[id]) * (current_time <= self$time), id))[id]
      invisible(self)
    },
    tilt_eps = function(trt, evnt, risk_evnt) {
      use <- data.frame(evnt = evnt,
                        offset = self$LH,
                        trt = trt,
                        Z1 = self$Z1,
                        Z0 = self$Z0,
                        risk_evnt = risk_evnt)
      self$eps <- check_na_coef(
        coef(
          speedglm::speedglm(
            evnt ~ 0 + offset(qlogis(offset)) + I(trt * Z1) + I((1 - trt) * Z0),
            family = binomial(),
            subset = risk_evnt == 1,
            data   = use
          )
        )
      )
      self$nuis$hzrd_on  <- bound01(plogis(qlogis(self$nuis$hzrd_on) + self$eps[1] * self$Z1))
      self$nuis$hzrd_off <- bound01(plogis(qlogis(self$nuis$hzrd_off) + self$eps[2] * self$Z0))
      invisible(self)
    },
    tilt_gamma = function(cens, risk_cens, all_time, trt) {
      use <- data.frame(cens = cens,
                        offset = self$GR,
                        H = self$H,
                        all_time = as.factor(all_time),
                        trt = trt,
                        risk_cens = risk_cens)
      self$gamma <- check_na_coef(
        coef(
          sw(speedglm::speedglm(
            cens ~ 0 + offset(qlogis(offset)) + H + all_time*trt,
            family = binomial(),
            subset = risk_cens == 1,
            data = use
          ))
        )
      )
      self$nuis$cens_on  <- bound01(plogis(qlogis(self$nuis$cens_on) + self$gamma[1] * self$H1))
      self$nuis$cens_off <- bound01(plogis(qlogis(self$nuis$cens_off) - self$gamma[1] * self$H0))
      invisible(self)
    },
    tilt_nu = function(trt, all_time) {
      use <- data.frame(trt = trt,
                        offset = self$nuis$trt_on,
                        M = self$M,
                        time = all_time)
      self$nu <- check_na_coef(
        coef(
          speedglm::speedglm(
            trt ~ 0 + offset(qlogis(offset)) + M,
            family = binomial(),
            subset = time == 1,
            data = use
          )
        )
      )
      self$nuis$trt_on  <- bound01(plogis(qlogis(self$nuis$trt_on) + self$nu * self$M))
      self$nuis$trt_off <- 1 - self$nuis$trt_on
      invisible(self)
    },
    update_crit = function(nobs) {
      self$iter <- self$iter + 1
      self$crit <- any(abs(c(self$eps, self$gamma[1], self$nu)) > 1e-3/nobs^(0.6))
      invisible(self)
    }
  )
)
