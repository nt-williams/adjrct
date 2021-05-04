Survival <- R6::R6Class(
  "Survival",
  public = list(
    outcome.formula = NULL,
    trt.formula = NULL,
    data = NULL,
    surv_data = NULL,
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
    initialize = function(outcome.formula, trt.formula, data) {
      self$outcome.formula <- outcome.formula
      self$trt.formula <- trt.formula
      self$data <- data
    },
    prepare_data = function(coarsen = 1) {
      self$trt <- get_time(self$trt.formula)

      if (!all(unique(self$data[[self$trt]]) %in% c(0, 1))) {
        stop("treatment should be coded as 0 and 1.", call. = FALSE)
      }

      self$time <- get_time(self$outcome.formula)
      self$status <- get_status(self$outcome.formula)
      self$covar <- get_covar(self$outcome.formula, self$trt)

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
    fit_nuis = function(algo, crossfit) {
      if (length(self$covar) < 2) algo <- "glm"

      on <- self$turn_on()
      off <- self$turn_off()

      fit_A <- glm(self$trt.formula, data = self$at_risk_trt(), family = "binomial")

      H_m <- model.matrix(self$formula_hzrd(), self$at_risk_evnt())
      R_m <- model.matrix(self$formula_cens(), self$at_risk_cens())

      H_mon <- model.matrix(self$formula_hzrd(), on)
      H_moff <- model.matrix(self$formula_hzrd(), off)

      R_mon <- model.matrix(self$formula_cens(), on)
      R_moff <- model.matrix(self$formula_cens(), off)

      if (algo == "glm") {
        fit_H <- fastglm::fastglm(x = H_m, y = as.matrix(self$at_risk_evnt()[["evnt"]]),
                                  family = "binomial")
        fit_R <- fastglm::fastglm(x = R_m, y = as.matrix(self$at_risk_cens()[["cens"]]),
                                  family = "binomial")

        self$nuisance <- list(
          hzrd_fit = fit_H,
          cens_fit = fit_R,
          trt_fit = fit_A,
          hzrd_off = bound01(predict(fit_H, newdata = H_moff, type = "response")),
          hzrd_on = bound01(predict(fit_H, newdata = H_mon, type = "response")),
          cens_off = bound01(predict(fit_R, newdata = R_moff, type = "response")),
          cens_on = bound01(predict(fit_R, newdata = R_mon, type = "response")),
          trt_off = bound01(1 - predict(fit_A, type = "response")),
          trt_on = bound01(predict(fit_A, type = "response"))
        )
        return(invisible(self))
      }

      if (algo == "lasso") {
        pen <- rep(1, ncol(R_m))
        pen[grep("^as.factor\\(all_time\\)", colnames(R_m))] <- 0

        folds <-
          foldsids(nrow(self$surv_data), self$surv_data[["survrctId"]],
                   automate_folds(nrow(self$surv_data)))
        fit_H <-
          glmnet::cv.glmnet(H_m,
                            as.matrix(self$at_risk_evnt()[["evnt"]]),
                            family = "binomial",
                            foldid = folds[self$risk_evnt == 1])
        fit_R <-
          glmnet::cv.glmnet(
            R_m,
            self$at_risk_cens()[["cens"]],
            family = "binomial",
            penalty.factor = pen,
            foldid = folds[self$risk_cens == 1]
          )

        self$nuisance <- list(
          hzrd_fit = fit_H,
          cens_fit = fit_R,
          trt_fit = fit_A,
          hzrd_off = bound01(as.vector(predict(fit_H, newx = H_moff, type = "response", s = "lambda.min"))),
          hzrd_on = bound01(as.vector(predict(fit_H, newx = H_mon, type = "response", s = "lambda.min"))),
          cens_off = bound01(as.vector(predict(fit_R, newx = R_moff, type = "response", s = "lambda.min"))),
          cens_on = bound01(as.vector(predict(fit_R, newx = R_mon, type = "response", s = "lambda.min"))),
          trt_off = bound01(as.vector(1 - predict(fit_A, type = "response"))),
          trt_on = bound01(as.vector(predict(fit_A, type = "response")))
        )
        return(invisible(self))
      }

      if (algo %in% c("rf", "xgboost")) {
        if (crossfit) {
          V <- automate_folds(nrow(self$surv_data))
          folds <- origami::make_folds(self$surv_data, cluster_ids = self$surv_data$survrctId, V = V)
          cf <- origami::cross_validate(cv_fun = cf_da_surv, algo = algo, folds = folds, self = self)
        } else {
          folds <- origami::make_folds(self$surv_data, cluster_ids = self$surv_data$survrctId, V = 1)
          folds[[1]]$training_set <- folds[[1]]$validation_set
          cf <- origami::cross_validate(cv_fun = cf_da_surv, algo = algo, folds = folds, self = self)
        }

        self$nuisance <- list(
          hzrd_fit = cf$hzrd_fit,
          cens_fit = cf$cens_fit,
          trt_fit = fit_A,
          hzrd_off = cf$hzrd_off[order(do.call("c", lapply(folds, function(x) x$validation_set)))],
          hzrd_on = cf$hzrd_on[order(do.call("c", lapply(folds, function(x) x$validation_set)))],
          cens_off = cf$cens_off[order(do.call("c", lapply(folds, function(x) x$validation_set)))],
          cens_on = cf$cens_on[order(do.call("c", lapply(folds, function(x) x$validation_set)))],
          trt_off = bound01(as.vector(1 - predict(fit_A, type = "response"))),
          trt_on = bound01(as.vector(predict(fit_A, type = "response")))
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
    formula_cens = function() {
      formula(paste("cens ~", self$trt, "* (", paste(c("as.factor(all_time)", self$covar), collapse = "+"), ")"))
    },
    formula_hzrd = function() {
      formula(paste("evnt ~", self$trt, "* (", paste(c("all_time", self$covar), collapse = "+"), ")"))
    },
    print = function(...) {
      cli::cli_text("{.strong survrct} metadata")
      cat("\n")
      cat("Outcome regression: ")
      print(self$outcome.formula)
      cat("        Propensity: ")
      print(self$trt.formula)
      cat("\n")
      cli::cli_ul(c("Estimate RMST with `rmst()`",
                    "Estimate survival probability with `survprob()`",
                    "Inspect nuisance parameter models with `get_fits()`"))
      cat("\n")
      cli::cli_text(cat("         "), "Estimator: TMLE")
      cli::cli_text(cat("   "), "Target variable: {self$trt}")
      cli::cli_text(cat("  "), "Status Indicator: {self$status}")
      cli::cli_text("Max coarsened time: {self$max_time}")
    }
  )
)

cf_da_surv <- function(algo, fold, self) {
  arH <- origami::training(self$risk_evnt)
  arR <- origami::training(self$risk_cens)
  train_H <- origami::training(self$surv_data)[arH == 1, ]
  train_R <- origami::training(self$surv_data)[arR == 1, ]
  folds_H <- origami::make_folds(train_H, cluster_ids = train_H[["survrctId"]], V = automate_folds(nrow(train_H)))
  folds_R <- origami::make_folds(train_R, cluster_ids = train_R[["survrctId"]], V = automate_folds(nrow(train_R)))
  on <- origami::validation(self$turn_on())
  off <- origami::validation(self$turn_off())

  method <- ifelse(algo == "rf", "ranger", "xgbTree")
  f_H <- reformulate(c("-1", self$trt, "all_time", self$covar))
  f_R <- reformulate(c("-1", self$trt, "as.factor(all_time)", self$covar))

  fit_H <- caret::train(
    x = model.matrix(f_H, data = train_H),
    y = factor(make.names(train_H[["evnt"]])),
    method = method,
    tuneLength = 5,
    trControl = caret::trainControl(
      method = "cv",
      classProbs = TRUE,
      search = "random",
      index = lapply(folds_H, function(x) x$training_set),
      indexOut = lapply(folds_H, function(x) x$validation_set)
    )
  )

  fit_R <- caret::train(
    x = model.matrix(f_R, data = train_R),
    y = factor(make.names(train_R[["cens"]])),
    method = method,
    tuneLength = 5,
    trControl = caret::trainControl(
      method = "cv",
      classProbs = TRUE,
      search = "random",
      index = lapply(folds_R, function(x) x$training_set),
      indexOut = lapply(folds_R, function(x) x$validation_set)
    )
  )

  list(
    hzrd_off = bound01(predict(fit_H, model.matrix(f_H, data = off), type = "prob")[, "X1"]),
    hzrd_on = bound01(predict(fit_H, model.matrix(f_H, data = on), type = "prob")[, "X1"]),
    cens_off = bound01(predict(fit_R, model.matrix(f_R, data = off), type = "prob")[, "X1"]),
    cens_on = bound01(predict(fit_R, model.matrix(f_R, data = on), type = "prob")[, "X1"]),
    hzrd_fit = fit_H,
    cens_fit = fit_R
  )
}
