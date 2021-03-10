Ordinal <- R6::R6Class(
  "Ordinal",
  public = list(
    formula = NULL,
    data = NULL,
    ordinal_data = NULL,
    estimator = NULL,
    trt = NULL,
    Y = NULL,
    covar = NULL,
    id = NULL,
    nobs = NULL,
    K = NULL,
    R = NULL,
    nuisance = list(),
    initialize = function(formula, target, data, estimator) {
      if (!all(unique(data[[target]]) %in% c(0, 1))) {
        stop("`target` should be coded as 0 and 1.", call. = FALSE)
      }

      if (!is.ordered(data[[all.vars(formula[[2]])]])) {
        stop("Outcome should be an ordered factor.", call. = FALSE)
      }
      self$formula <- formula
      self$data <- as.data.frame(data)
      self$trt <- target
      self$estimator <- estimator
    },
    prepare_data = function() {
      self$Y <- all.vars(self$formula[[2]])
      self$covar <- get_covar(self$formula, self$trt)
      self$K <- length(unique(self$data[[self$Y]]))
      self$nobs <- nrow(self$data)
      self$id <- rep(1:self$nobs, each = self$K - 1)

      W <- self$data[self$id, self$covar, drop = FALSE]
      A <- self$data[[self$trt]][self$id]
      kl <- as.factor(rep(1:(self$K - 1), self$nobs))
      Yl <- as.numeric(as.numeric(self$data[[self$Y]][self$id]) == kl)
      self$R <- unlist(lapply(tapply(Yl, self$id, cumsum), function(x) as.numeric(cumsum(x) <= 1)))

      self$ordinal_data <- data.frame(ordinalrctId = self$id, W, A, kl, Y = Yl)
      invisible(self)
    },
    fit_nuis = function(lasso) {
      if (length(self$covar) < 2) lasso <- FALSE
      if (!lasso) {
        fit_Q <- glm(self$formula_y(), family = binomial(), data = self$at_risk())
        fit_A <- glm(self$formula_trt(), family = binomial(), data = self$data)

        self$nuisance <- list(
          hzrd_fit = fit_Q,
          trt_fit = fit_A,
          H_off = as.vector(predict(fit_Q, newdata = self$turn_off(), type = "response")),
          H_on = as.vector(predict(fit_Q, newdata = self$turn_on(), type = "response")),
          trt_off = as.vector(1 - predict(fit_A, type = "response")),
          trt_on = as.vector(predict(fit_A, type = "response"))
        )
        return(invisible(self))
      }
      on <- self$turn_on()
      off <- self$turn_off()

      H_m <- model.matrix(self$formula_y(), self$at_risk())[, -1, drop = FALSE]
      A_m <- model.matrix(self$formula_trt(), self$data)[, -1, drop = FALSE]

      H_mon <- model.matrix(self$formula_y(), on)[, -1, drop = FALSE]
      H_moff <- model.matrix(self$formula_y(), off)[, -1, drop = FALSE]

      A_o <- model.matrix(self$formula_trt(), self$data)[, -1, drop = FALSE]

      folds <- foldsids(nrow(self$ordinal_data), self$ordinal_data[["ordinalrctId"]], 10)

      fit_Q <- glmnet::cv.glmnet(H_m, as.matrix(self$at_risk()[["Y"]]),
                                 family = "binomial", foldid = folds[self$R == 1])
      fit_A <- glmnet::cv.glmnet(A_m, self$data[[self$trt]],
                                 family = "binomial", foldid = as.vector(tapply(folds, self$id, function(x) x[1])))

      self$nuisance <- list(
        hzrd_fit = fit_Q,
        trt_fit = fit_A,
        H_off = bound01(as.vector(predict(fit_Q, newx = H_moff, type = "response", s = "lambda.min"))),
        H_on = bound01(as.vector(predict(fit_Q, newx = H_mon, type = "response", s = "lambda.min"))),
        trt_off = bound01(as.vector(1 - predict(fit_A, newx = A_o, type = "response", s = "lambda.min"))),
        trt_on = bound01(as.vector(predict(fit_A, newx = A_o, type = "response", s = "lambda.min")))
      )
      invisible(self)
    },
    at_risk = function() {
      self$ordinal_data[self$R == 1, ]
    },
    turn_on = function() {
      out <- self$ordinal_data
      out[[self$trt]] <- rep(1, nrow(self$ordinal_data))
      return(out)
    },
    turn_off = function() {
      out <- self$ordinal_data
      out[[self$trt]] <- rep(0, nrow(self$ordinal_data))
      return(out)
    },
    formula_trt = function() {
      if (length(self$covar) == 0) {
        covar <- 1
      } else {
        covar <- self$covar
      }
      formula(paste(self$trt, "~", paste(covar, collapse = "+")))
    },
    formula_y = function() {
      formula(paste("Y ~ -1 + kl*(", paste(c("A", self$covar), collapse = "+"), ")"))
    },
    print = function(...) {
      cli::cli_text("{.strong ordinalrct} metadata")
      cat("\n")
      print(self$formula)
      cat("\n")
      cli::cli_ul(c("Estimate log odds ratio with `log_or()`",
                    "Estimate Mann-Whitney with `mannwhitney()`",
                    "Estimate with `cdf()`",
                    "Estimate with `pmf()`",
                    "Inspect nuisance parameter models with `get_fits()`"))
      cat("\n")
      cli::cli_text(cat("         "), "Estimator: {self$estimator}")
      cli::cli_text(cat("   "), "Target variable: {self$trt}")
      cli::cli_text(cat("  "), "Outcome variable: {self$Y}")
      cli::cli_text(cat("    "), "Adjustment set: {self$covar}")
    }
  )
)
