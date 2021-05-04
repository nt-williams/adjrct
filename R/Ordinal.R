Ordinal <- R6::R6Class(
  "Ordinal",
  public = list(
    outcome.formula = NULL,
    trt.formula = NULL,
    data = NULL,
    ordinal_data = NULL,
    trt = NULL,
    Y = NULL,
    covar = NULL,
    id = NULL,
    nobs = NULL,
    K = NULL,
    R = NULL,
    nuisance = list(),
    initialize = function(outcome.formula, trt.formula, data) {
      self$outcome.formula <- outcome.formula
      self$data <- as.data.frame(data)
      self$trt.formula <- trt.formula
    },
    prepare_data = function() {
      self$Y <- get_time(self$outcome.formula)
      self$trt <- get_time(self$trt.formula)

      if (!all(unique(self$data[[self$trt]]) %in% c(0, 1))) {
        stop("The treatment should be coded as 0 and 1.", call. = FALSE)
      }

      if (!is.ordered(self$data[[self$Y]])) {
        stop("The outcome should be an ordered factor.", call. = FALSE)
      }

      self$covar <- get_covar(self$outcome.formula, self$trt)
      self$K <- length(unique(self$data[[self$Y]]))
      self$nobs <- nrow(self$data)
      self$id <- rep(1:self$nobs, each = self$K - 1)
      W <- self$data[self$id, self$covar, drop = FALSE]
      A <- self$data[[self$trt]][self$id]
      kl <- as.factor(rep(1:(self$K - 1), self$nobs))
      Yl <- as.numeric(as.numeric(self$data[[self$Y]][self$id]) == kl)
      self$R <- unlist(lapply(tapply(Yl, self$id, cumsum), function(x) as.numeric(cumsum(x) <= 1)))

      self$ordinal_data <- data.frame(ordinalrctId = self$id, atRisk = self$R, W, A, kl, Y = Yl)
      invisible(self)
    },
    fit_nuis = function(algo, crossfit = TRUE) {
      if (length(self$covar) < 2) algo <- "glm"

      fit_A <- glm(self$trt.formula, family = binomial(), data = self$data)

      if (algo == "glm") {
        fit_Q <- glm(self$formula_y(), family = binomial(), data = self$at_risk())

        self$nuisance <- list(
          hzrd_fit = fit_Q,
          trt_fit = fit_A,
          H_off = bound01(as.vector(predict(fit_Q, newdata = self$turn_off(), type = "response"))),
          H_on = bound01(as.vector(predict(fit_Q, newdata = self$turn_on(), type = "response"))),
          trt_off = bound01(as.vector(1 - predict(fit_A, type = "response"))),
          trt_on = bound01(as.vector(predict(fit_A, type = "response")))
        )
        return(invisible(self))
      }

      on <- self$turn_on()
      off <- self$turn_off()

      if (algo == "lasso") {
        H_m <- model.matrix(self$formula_y(), self$at_risk())
        H_mon <- model.matrix(self$formula_y(), on)
        H_moff <- model.matrix(self$formula_y(), off)

        folds <- foldsids(nrow(self$ordinal_data), self$ordinal_data[["ordinalrctId"]], 10)
        fit_Q <- glmnet::cv.glmnet(H_m, as.matrix(self$at_risk()[["Y"]]),
                                   family = "binomial", foldid = folds[self$R == 1])

        self$nuisance <- list(
          hzrd_fit = fit_Q,
          trt_fit = fit_A,
          H_off = bound01(as.vector(predict(fit_Q, newx = H_moff, type = "response", s = "lambda.min"))),
          H_on = bound01(as.vector(predict(fit_Q, newx = H_mon, type = "response", s = "lambda.min"))),
          trt_off = bound01(as.vector(1 - predict(fit_A, type = "response"))),
          trt_on = bound01(as.vector(predict(fit_A, type = "response")))
        )
        return(invisible(self))
      }

      if (algo %in% c("rf", "xgboost", "earth")) {
        if (crossfit) {
          V <- automate_folds(nrow(self$ordinal_data))
          folds <- origami::make_folds(self$ordinal_data, cluster_ids = self$ordinal_data$ordinalrctId, V = V)
          cf <- origami::cross_validate(cv_fun = cf_da_ord, algo = algo, folds = folds, self = self)
        } else {
          folds <- origami::make_folds(self$ordinal_data, cluster_ids = self$ordinal_data$ordinalrctId, V = 1)
          folds[[1]]$training_set <- folds[[1]]$validation_set
          cf <- origami::cross_validate(cv_fun = cf_da_ord, algo = algo, folds = folds, self = self)
        }

        self$nuisance <- list(
          hzrd_fit = cf$hzrd_fit,
          trt_fit = fit_A,
          H_off = cf$H_off[order(do.call("c", lapply(folds, function(x) x$validation_set)))],
          H_on = cf$H_on[order(do.call("c", lapply(folds, function(x) x$validation_set)))],
          trt_off = bound01(as.vector(1 - predict(fit_A, type = "response"))),
          trt_on = bound01(as.vector(predict(fit_A, type = "response")))
        )
        return(invisible(self))
      }
    },
    at_risk = function() {
      self$ordinal_data[self$R == 1, ]
    },
    turn_on = function() {
      out <- self$ordinal_data
      out[["A"]] <- rep(1, nrow(self$ordinal_data))
      return(out)
    },
    turn_off = function() {
      out <- self$ordinal_data
      out[["A"]] <- rep(0, nrow(self$ordinal_data))
      return(out)
    },
    formula_y = function() {
      formula(paste("Y ~ -1 + kl*(", paste(c("A", self$covar), collapse = "+"), ")"))
    },
    print = function(...) {
      cli::cli_text("{.strong ordinalrct} metadata")
      cat("\n")
      cat("Outcome regression: ")
      print(self$outcome.formula)
      cat("        Propensity: ")
      print(self$trt.formula)
      cat("\n")
      cli::cli_ul(c("Estimate log odds ratio with `log_or()`",
                    "Estimate Mann-Whitney with `mannwhitney()`",
                    "Estimate with `cdf()`",
                    "Estimate with `pmf()`",
                    "Inspect nuisance parameter models with `get_fits()`"))
      cat("\n")
      cli::cli_text(cat("         "), "Estimator: TMLE")
      cli::cli_text(cat("   "), "Target variable: {self$trt}")
      cli::cli_text(cat("  "), "Outcome variable: {self$Y}")
    }
  )
)

cf_da_ord <- function(algo, fold, self) {
  train <- origami::training(self$ordinal_data)
  train <- train[train$atRisk == 1, ]
  folds <- origami::make_folds(train, cluster_ids = train[["ordinalrctId"]], V = automate_folds(nrow(train)))
  f <- reformulate(c("-1", "kl", self$trt, self$covar))

  fit <- caret::train(
    x = model.matrix(f, data = train),
    y = factor(make.names(train[["Y"]])),
    method = method <- ifelse(algo == "rf", "ranger", "xgbTree"),
    tuneLength = 5,
    trControl = caret::trainControl(
      method = "cv",
      classProbs = TRUE,
      search = "random",
      index = lapply(folds, function(x) x$training_set),
      indexOut = lapply(folds, function(x) x$validation_set)
    )
  )

  list(H_off = bound01(predict(fit, model.matrix(f, data = origami::validation(self$turn_off())),
                               type = "prob")[, "X1"]),
       H_on = bound01(predict(fit, model.matrix(f, data = origami::validation(self$turn_on())),
                              type = "prob")[, "X1"]),
       hzrd_fit = fit)
}
