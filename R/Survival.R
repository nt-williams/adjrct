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
    fit_nuis = function(algo, crossfit) {
      if (length(self$covar) < 2) algo <- "glm"

      if (algo == "glm") {
        fit_H <- glm(self$formula_hzrd(), data = self$at_risk_evnt(), family = "binomial")
        fit_R <- glm(self$formula_cens(), data = self$at_risk_cens(), family = "binomial")
        fit_A <- glm(self$formula_trt(), data = self$at_risk_trt(), family = "binomial")

        self$nuisance <- list(
          hzrd_fit = fit_H,
          cens_fit = fit_R,
          trt_fit = fit_A,
          hzrd_off = bound01(predict(fit_H, newdata = self$turn_off(), type = "response")),
          hzrd_on = bound01(predict(fit_H, newdata = self$turn_on(), type = "response")),
          cens_off = bound01(predict(fit_R, newdata = self$turn_off(), type = "response")),
          cens_on = bound01(predict(fit_R, newdata = self$turn_on(), type = "response")),
          trt_off = bound01(1 - predict(fit_A, type = "response")),
          trt_on = bound01(predict(fit_A, type = "response"))
        )
        return(invisible(self))
      }

      on <- self$turn_on()
      off <- self$turn_off()

      if (algo == "lasso") {
        H_m <- model.matrix(self$formula_hzrd(), self$at_risk_evnt())[, -1]
        H_mon <- model.matrix(self$formula_hzrd(), on)[, -1]
        H_moff <- model.matrix(self$formula_hzrd(), off)[, -1]
        folds <- foldsids(nrow(self$surv_data), self$surv_data[["survrctId"]], 10)
        fit_H <- glmnet::cv.glmnet(H_m, as.matrix(self$at_risk_evnt()[["evnt"]]),
                                   family = "binomial", foldid = folds[self$risk_evnt == 1])

        # for the paper, only the hazard will be estimated using variable selection
        fit_R <- glm(self$formula_cens(), data = self$at_risk_cens(), family = binomial())
        fit_A <- glm(self$formula_trt(), family = binomial(), data = self$at_risk_trt())
        self$nuisance <- list(
          hzrd_fit = fit_H,
          cens_fit = fit_R,
          trt_fit = fit_A,
          hzrd_off = bound01(as.vector(predict(fit_H, newx = H_moff, type = "response", s = "lambda.min"))),
          hzrd_on = bound01(as.vector(predict(fit_H, newx = H_mon, type = "response", s = "lambda.min"))),
          cens_off = bound01(predict(fit_R, newdata = self$turn_off(), type = "response")),
          cens_on = bound01(predict(fit_R, newdata = self$turn_on(), type = "response")),
          trt_off = bound01(as.vector(1 - predict(fit_A, type = "response"))),
          trt_on = bound01(as.vector(predict(fit_A, type = "response")))
        )
        return(invisible(self))
      }

      if (algo %in% c("rf", "xgboost", "earth")) {
        if (crossfit) {
          V <- automate_folds(nrow(self$surv_data))
          folds <- origami::make_folds(self$surv_data, cluster_ids = self$surv_data$survrctId, V = V)
          cf <- origami::cross_validate(cv_fun = cf_da_surv, algo = algo, folds = folds, self = self)
        } else {
          folds <- origami::make_folds(self$surv_data, cluster_ids = self$surv_data$survrctId, V = 1)
          folds[[1]]$training_set <- folds[[1]]$validation_set
          cf <- origami::cross_validate(cv_fun = cf_da_surv, algo = algo, folds = folds, self = self)
        }

        # for the paper, only the hazard will be estimated using variable selection
        fit_R <- glm(self$formula_cens(), data = self$at_risk_cens(), family = binomial())
        fit_A <- glm(self$formula_trt(), family = binomial(), data = self$at_risk_trt())
        self$nuisance <- list(
          hzrd_fit = cf$hzrd_fit,
          cens_fit = fit_R,
          trt_fit = fit_A,
          hzrd_off = cf$hzrd_off[order(do.call("c", lapply(folds, function(x) x$validation_set)))],
          hzrd_on = cf$hzrd_on[order(do.call("c", lapply(folds, function(x) x$validation_set)))],
          cens_off = bound01(predict(fit_R, newdata = self$turn_off(), type = "response")),
          cens_on = bound01(predict(fit_R, newdata = self$turn_on(), type = "response")),
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
    formula_trt = function() {
      formula(paste(self$trt, "~", 1))
    },
    formula_cens = function() {
      formula(paste("cens ~", self$trt, "* as.factor(all_time)"))
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

cf_da_surv <- function(algo, fold, self) {
  arH <- origami::training(self$risk_evnt)
  train_H <- origami::training(self$surv_data)[arH == 1, ]
  VH <- automate_folds(nrow(train_H))

  on <- origami::validation(self$turn_on())
  off <- origami::validation(self$turn_off())

  if (algo == "rf") {
    fit_H <- caret::train(
      x = train_H[, c("all_time", self$trt, self$covar)],
      y = factor(make.names(train_H[["evnt"]])),
      method = "ranger",
      tuneLength = 5,
      trControl = caret::trainControl(
        method = "cv",
        classProbs = TRUE,
        search = "random",
        index = lapply(1:VH, function(x) which(x != foldsids(nrow(train_H), train_H[["survrctId"]], VH)))
      )
    )

    out <- list(hzrd_off = bound01(predict(fit_H, off, type = "prob")[, "X1"]),
                hzrd_on = bound01(predict(fit_H, on, type = "prob")[, "X1"]),
                hzrd_fit = fit_H)
    return(out)
  }

  if (algo == "xgboost") {
    fit_H <- caret::train(
      x = model.matrix(reformulate(c("all_time", self$trt, self$covar)), data = train_H),
      y = factor(make.names(train_H[["evnt"]])),
      method = "xgbTree",
      tuneLength = 5,
      trControl = caret::trainControl(
        method = "cv",
        classProbs = TRUE,
        search = "random",
        index = lapply(1:VH, function(x) which(x != foldsids(nrow(train_H), train_H[["survrctId"]], VH)))
      )
    )

    out <- list(hzrd_off = bound01(predict(fit_H, model.matrix(reformulate(c("all_time", self$trt, self$covar)), data = off), type = "prob")[, "X1"]),
                hzrd_on = bound01(predict(fit_H, model.matrix(reformulate(c("all_time", self$trt, self$covar)), data = on), type = "prob")[, "X1"]),
                hzrd_fit = fit_H)
    return(out)
  }

  if (algo == "earth") {
    fit_H <- caret::train(
      x = train_H[, c("all_time", self$trt, self$covar)],
      y = factor(make.names(train_H[["evnt"]])),
      method = "earth",
      tuneLength = 5,
      trControl = caret::trainControl(
        method = "cv",
        classProbs = TRUE,
        search = "random",
        index = lapply(1:VH, function(x) which(x != foldsids(nrow(train_H), train_H[["survrctId"]], VH)))
      )
    )

    out <- list(hzrd_off = bound01(predict(fit_H, off, type = "prob")[, "X1"]),
                hzrd_on = bound01(predict(fit_H, on, type = "prob")[, "X1"]),
                hzrd_fit = fit_H)
    return(out)
  }
}
