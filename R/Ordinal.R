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
      self$data <- data
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

      self$ordinal_data <- data.frame(W, A, kl, Y = Yl)
      invisible(self)
    },
    fit_nuis = function() {
      fit_Q <- glm(self$formula_y(), family = binomial(), data = self$at_risk())
      fit_A <- glm(self$formula_trt(), family = binomial(), data = self$data)

      self$nuisance <- list(
        fit_H = fit_Q,
        fit_A = fit_A,
        H_off = as.vector(predict(fit_Q, newdata = self$turn_off(), type = "response")),
        H_on = as.vector(predict(fit_Q, newdata = self$turn_on(), type = "response")),
        trt_off = as.vector(1 - predict(fit_A, type = "response")),
        trt_on = as.vector(predict(fit_A, type = "response"))
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
      if (is.null(self$covar)) {
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
                    "Get CDF with `cdf()`",
                    "Inspect nuisance parameter models with `get_fits()`"))
      cat("\n")
      cli::cli_text(cat("         "), "Estimator: {self$estimator}")
      cli::cli_text(cat("   "), "Target variable: {self$trt}")
      cli::cli_text(cat("  "), "Outcome variable: {self$Y}")
      cli::cli_text(cat("    "), "Adjustment set: {self$covar}")
    }
  )
)
