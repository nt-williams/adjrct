Survival <- R6::R6Class(
  "Survival",
  public = list(
    formula = NULL,
    data = NULL,
    surv_data = NULL,
    estimator = NULL,
    trt = NULL,
    status = NULL,
    covar = NULL,
    time = NULL,
    horizon = list(),
    id = NULL,
    nobs = NULL,
    all_time = NULL,
    max_time = NULL,
    risk_evnt = NULL,
    risk_cens = NULL,
    nuisance = list(),
    initialize = function(formula, target, data, estimator) {
      if (!all(unique(data[[target]]) %in% c(0, 1))) {
        stop("`target` should be coded as 0 and 1.", call. = FALSE)
      }
      self$formula <- formula
      self$data <- data
      self$trt <- target
      self$estimator <- estimator
    },
    prepare_data = function(coarsen = 1) {
      # TODO check_status() i.e., status is 0 and 1

      self$time <- get_time(self$formula)
      self$status <- get_status(self$formula)
      self$covar <- get_covar(self$formula, self$trt)

      if (coarsen > 1) self$data[[self$time]] <- self$data[[self$time]] %/% coarsen + 1

      self$nobs <- nrow(self$data)
      self$max_time <- max(self$data[[self$time]])
      self$all_time <- rep(1:self$max_time, self$nobs)
      evnt <- cens <- rep(NA, self$nobs*self$max_time)
      self$risk_evnt <- self$risk_cens <- 1*(self$all_time == 1)

      for (t in 1:self$max_time) {
        cens[self$all_time == t] <- (1 - self$data[[self$status]]) * (self$data[[self$time]] == t)
        evnt[self$all_time == t] <- self$data[[self$status]] * (self$data[[self$time]] == t)
        self$risk_evnt[self$all_time == t] <- (self$data[[self$time]] >= t)
        self$risk_cens[self$all_time == t] <- (self$data[[self$time]] > t) * self$data[[self$status]] +
          (self$data[[self$time]] >= t) * (1 - self$data[[self$status]])
      }

      self$surv_data <- data.frame(survrctId = rep(1:self$nobs, each = self$max_time),
                                   self$data[as.numeric(gl(self$nobs, self$max_time)), ],
                                   all_time = self$all_time, evnt, cens)
      invisible(self)
    },
    fit_nuis = function(lasso) {
      self$nuisance <- nuisance(self, lasso)
      invisible(self)
    },
    evaluate_horizon = function(horizon = NULL, estimand) {
      # also need to add check for the minimum time
      if (!is.null(horizon)) {
        if (max(horizon) > self$max_time) {
          stop("Horizon is greater than max observed time!", call. = FALSE)
        }
        self$horizon <- horizon
      } else {
        if (estimand == "rmst") {
          self$horizon <- 2:self$max_time
        } else {
          self$horizon <- 1:self$max_time
        }
      }
      invisible(self)
    },
    at_risk_evnt = function() {
      self$surv_data[self$risk_evnt == 1, ]
    },
    at_risk_cens = function() {
      self$surv_data[self$risk_cens == 1, ]
    },
    at_risk_trt = function() {
      self$surv_data[self$all_time == 1, ]
    },
    turn_on = function() {
      out <- self$surv_data
      out[[self$trt]] <- rep(1, nrow(self$surv_data))
      return(out)
    },
    turn_off = function() {
      out <- self$surv_data
      out[[self$trt]] <- rep(0, nrow(self$surv_data))
      return(out)
    },
    predictors = function() {
      c("all_time", self$trt, self$covar)
    },
    get_var = function(var) {
      self$surv_data[[self[[var]]]]
    },
    time_indicator = function() {
      outer(self$all_time, 1:self$max_time, "<=")
    },
    formula_trt = function() {
      if (length(self$covar) == 0) {
        covar <- 1
      } else {
        covar <- self$covar
      }
      formula(paste(self$trt, "~", paste(covar, collapse = "+")))
    },
    formula_cens = function() {
      if (self$estimator %in% c("tmle", "aipw")) {
        formula(paste("cens ~", self$trt, "* (", paste(c("as.factor(all_time)", self$covar), collapse = "+"), ")"))
      } else {
        formula(paste("cens ~ as.factor(all_time) * ", self$trt))
      }
    },
    formula_hzrd = function() {
      if (self$estimator %in% c("tmle", "aipw")) {
        formula(paste("evnt ~", self$trt, "* (", paste(c("all_time", self$covar), collapse = "+"), ")"))
      } else {
        formula(paste("evnt ~ all_time * ", self$trt))
      }
    },
    print = function(...) {
      cli::cli_text("{.strong survrct} metadata")
      cat("\n")
      print(self$formula)
      cat("\n")
      cli::cli_ul(c("Estimate RMST with `rmst()`",
                    "Estimate survival probability with `survprob()`",
                    "Inspect nuisance parameter models with `get_fits()`"))
      cat("\n")
      cli::cli_text(cat("         "), "Estimator: {self$estimator}")
      cli::cli_text(cat("   "), "Target variable: {self$trt}")
      cli::cli_text(cat("  "), "Status Indicator: {self$status}")
      cli::cli_text(cat("    "), "Adjustment set: {self$covar}")
      cli::cli_text("Max coarsened time: {self$max_time}")
    }
  )
)
