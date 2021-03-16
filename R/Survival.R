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
      self$time <- get_time(self$formula)
      self$status <- get_status(self$formula)
      self$covar <- get_covar(self$formula, self$trt)

      if (!all(unique(self$data[[self$status]]) %in% c(0, 1))) {
        stop("`status` should be coded as 0 and 1.", call. = FALSE)
      }

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
    fit_nuis = function(algo) {
      if (length(self$covar) < 2) algo <- "glm"

      if (algo == "glm") {
        fit_H <- glm(self$formula_hzrd(), data = self$at_risk_evnt(), family = "binomial")
        fit_R <- glm(self$formula_cens(), data = self$at_risk_cens(), family = "binomial")
        fit_A <- glm(self$formula_trt(), data = self$at_risk_trt(), family = "binomial")

        self$nuisance <- list(
          hzrd_fit = fit_H,
          cens_fit = fit_R,
          trt_fit = fit_A,
          hzrd_off = predict(fit_H, newdata = self$turn_off(), type = "response"),
          hzrd_on = predict(fit_H, newdata = self$turn_on(), type = "response"),
          cens_off = predict(fit_R, newdata = self$turn_off(), type = "response"),
          cens_on = predict(fit_R, newdata = self$turn_on(), type = "response"),
          trt_off = 1 - predict(fit_A, type = "response"),
          trt_on = predict(fit_A, type = "response")
        )
        return(invisible(self))
      }

      on <- self$turn_on()
      off <- self$turn_off()

      if (algo == "lasso") {
        H_m <- model.matrix(self$formula_hzrd(), self$at_risk_evnt())[, -1]
        R_m <- model.matrix(self$formula_cens(), self$at_risk_cens())[, -1]
        A_m <- model.matrix(self$formula_trt(), self$at_risk_trt())[, -1]
        H_mon <- model.matrix(self$formula_hzrd(), on)[, -1]
        H_moff <- model.matrix(self$formula_hzrd(), off)[, -1]
        R_mon <- model.matrix(self$formula_cens(), on)[, -1]
        R_moff <- model.matrix(self$formula_cens(), off)[, -1]
        pen <- rep(1, ncol(R_m))
        pen[grep("^as.factor\\(all_time\\)", colnames(R_m))] <- 0
        folds <- foldsids(nrow(self$surv_data), self$surv_data[["survrctId"]], 10)
        fit_H <- glmnet::cv.glmnet(H_m, as.matrix(self$at_risk_evnt()[["evnt"]]),
                                   family = "binomial", foldid = folds[self$risk_evnt == 1])
        fit_R <- glmnet::cv.glmnet(R_m, self$at_risk_cens()[["cens"]],
                                   family = "binomial", penalty.factor = pen,
                                   foldid = folds[self$risk_cens == 1])
        fit_A <- glmnet::cv.glmnet(A_m, self$at_risk_trt()[[self$trt]],
                                   family = "binomial", foldid = folds[self$all_time == 1])
        self$nuisance <- list(
          hzrd_fit = fit_H,
          cens_fit = fit_R,
          trt_fit = fit_A,
          hzrd_off = bound01(as.vector(predict(fit_H, newx = H_moff, type = "response", s = "lambda.min"))),
          hzrd_on = bound01(as.vector(predict(fit_H, newx = H_mon, type = "response", s = "lambda.min"))),
          cens_off = bound01(as.vector(predict(fit_R, newx = R_moff, type = "response", s = "lambda.min"))),
          cens_on = bound01(as.vector(predict(fit_R, newx = R_mon, type = "response", s = "lambda.min"))),
          trt_off = bound01(as.vector(1 - predict(fit_A, newx = A_m, type = "response", s = "lambda.min"))),
          trt_on = bound01(as.vector(predict(fit_A, newx = A_m, type = "response", s = "lambda.min")))
        )
        return(invisible(self))
      }

      if (algo == "rf") {
        fit_H <- caret::train(
          x = self$at_risk_evnt()[, c("all_time", self$trt, self$covar)],
          y = factor(make.names(self$at_risk_evnt()[["evnt"]])),
          method = "ranger",
          tuneLength = 5,
          trControl = caret::trainControl(
            method = "cv",
            classProbs = TRUE,
            search = "random",
            index = lapply(1:10, function(x) which(x != foldsids(nrow(self$at_risk_evnt()), self$at_risk_evnt()[["survrctId"]], 10)))
          )
        )
        fit_R <- caret::train(
          x = self$at_risk_cens()[, c("all_time", self$trt, self$covar)],
          y = factor(make.names(self$at_risk_cens()[["cens"]])),
          method = "ranger",
          tuneLength = 5,
          trControl = caret::trainControl(
            method = "cv",
            classProbs = TRUE,
            search = "random",
            index = lapply(1:10, function(x) which(x != foldsids(nrow(self$at_risk_cens()), self$at_risk_cens()[["survrctId"]], 10)))
          )
        )
        fit_A <- glm(formula(paste(self$trt, "~", 1)), family = binomial(), data = self$at_risk_trt())
        self$nuisance <- list(
          hzrd_fit = fit_H,
          cens_fit = fit_R,
          trt_fit = fit_A,
          hzrd_off = bound01(predict(fit_H, off, type = "prob")[, "X1"]),
          hzrd_on = bound01(predict(fit_H, on, type = "prob")[, "X1"]),
          cens_off = bound01(predict(fit_R, off, type = "prob")[, "X1"]),
          cens_on = bound01(predict(fit_R, on, type = "prob")[, "X1"]),
          trt_off = as.vector(1 - predict(fit_A, type = "response")),
          trt_on = as.vector(predict(fit_A, type = "response"))
        )
        return(invisible(self))
      }
    },
    evaluate_horizon = function(horizon = NULL, estimand) {
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
      formula(paste("cens ~", self$trt, "* (", paste(c("as.factor(all_time)", self$covar), collapse = "+"), ")"))
    },
    formula_hzrd = function() {
      formula(paste("evnt ~", self$trt, "* (", paste(c("all_time", self$covar), collapse = "+"), ")"))
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
