#' @export
print.rmst <- function(x, ...) {
  cli::cli_text("{.strong RMST Estimator}: {x$estimator}")
  if (length(x$horizon) > 1) {
    cat("\n")
    cat("           Confidence level: 95%\n")
    cat("          Mult. Bootstrap C:", x$estimates$mbcv_theta, "\n")
    cat("\n")
    print(head(format_est(x)))
    cli::cli_text(cli::col_red("Access all estimates with `all_estimates()`"))
  } else {
    cli::cli_text(cat("  "), "{.strong Time horizon}: {x$horizon}")
    cat("\n")
    cli::cli_text("{.strong Arm-specific RMST:}")
    cli::cli_text(cli::col_blue(cli::style_italic("Treatment Arm")))
    cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates[[1]]$arm1, 2)}")
    cli::cli_text(cat("    "), "{.strong Std. error}: {round(x$estimates[[1]]$arm1.std.error, 2)}")
    cli::cli_text(cat("        "), "{.strong 95% CI}: ({round(x$estimates[[1]]$arm1.conf.low, 2)}, {round(x$estimates[[1]]$arm1.conf.high, 2)})")
    cli::cli_text(cli::col_red(cli::style_italic("Control Arm")))
    cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates[[1]]$arm0, 2)}")
    cli::cli_text(cat("    "), "{.strong Std. error}: {round(x$estimates[[1]]$arm0.std.error, 2)}")
    cli::cli_text(cat("        "), "{.strong 95% CI}: ({round(x$estimates[[1]]$arm0.conf.low, 2)}, {round(x$estimates[[1]]$arm0.conf.high, 2)})")
    cat("\n")
    cli::cli_text("{.strong Treatment Effect:}")
    cli::cli_text(cli::col_green(cli::style_italic("Additive effect")))
    cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates[[1]]$theta, 2)}")
    cli::cli_text(cat("    "), "{.strong Std. error}: {round(x$estimates[[1]]$std.error, 2)}")
    cli::cli_text(cat("        "), "{.strong 95% CI}: ({round(x$estimates[[1]]$theta.conf.low, 2)}, {round(x$estimates[[1]]$theta.conf.high, 2)})")
  }
}

#' @export
print.survprob <- function(x, ...) {
  cli::cli_text("{.strong Survival Probability Estimator}: {x$estimator}")
  if (length(x$horizon) > 1) {
    cat("\n")
    cat("           Confidence level: 95%\n")
    cat("          Mult. Bootstrap C:", x$estimates$mbcv_theta, "\n")
    cat("\n")
    print(head(format_est(x)))
    cli::cli_text(cli::col_red("Access all estimates with `all_estimates()`"))
  } else {
    cli::cli_text(cat("  "), "{.strong Time horizon}: {x$horizon}")
    cat("\n")
    cli::cli_text("{.strong Arm-specific Survival Probability:}")
    cli::cli_text(cli::col_blue(cli::style_italic("Treatment Arm")))
    cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates[[1]]$arm1, 2)}")
    cli::cli_text(cat("    "), "{.strong Std. error}: {round(x$estimates[[1]]$arm1.std.error, 2)}")
    cli::cli_text(cat("        "), "{.strong 95% CI}: ({round(x$estimates[[1]]$arm1.conf.low, 2)}, {round(x$estimates[[1]]$arm1.conf.high, 2)})")
    cli::cli_text(cli::col_red(cli::style_italic("Control Arm")))
    cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates[[1]]$arm0, 2)}")
    cli::cli_text(cat("    "), "{.strong Std. error}: {round(x$estimates[[1]]$arm0.std.error, 2)}")
    cli::cli_text(cat("        "), "{.strong 95% CI}: ({round(x$estimates[[1]]$arm0.conf.low, 2)}, {round(x$estimates[[1]]$arm0.conf.high, 2)})")
    cat("\n")
    cli::cli_text("{.strong Treatment Effect:}")
    cli::cli_text(cli::col_green(cli::style_italic("Additive effect")))
    cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates[[1]]$theta, 2)}")
    cli::cli_text(cat("    "), "{.strong Std. error}: {round(x$estimates[[1]]$std.error, 2)}")
    cli::cli_text(cat("        "), "{.strong 95% CI}: ({round(x$estimates[[1]]$theta.conf.low, 2)}, {round(x$estimates[[1]]$theta.conf.high, 2)})")
  }
}

#' @export
print.lor <- function(x, ...) {
  cli::cli_text("{.strong Log OR Estimator}: {x$estimator}")
    cat("\n")
    cli::cli_text("{.strong Arm-specific log odds:}")
    cli::cli_text(cli::col_blue(cli::style_italic("Treatment Arm")))
    cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates$lor$arm1, 2)}")
    cli::cli_text(cli::col_red(cli::style_italic("Control Arm")))
    cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates$lor$arm0, 2)}")
    cat("\n")
    cli::cli_text("{.strong Average log odds ratio:}")
    cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates$lor$theta, 2)}")
    cli::cli_text(cat("    "), "{.strong Std. error}: {round(x$estimates$std.error, 2)}")
    cli::cli_text(cat("        "), "{.strong 95% CI}: ({round(x$estimates$ci[1], 2)}, {round(x$estimates$ci[2], 2)})")
}

#' @export
print.cdf <- function(x, ...) {
  cli::cli_text("{.strong CDF Estimator}: {x$estimator}")
  cat("\n")
  cli::cli_text("{.strong Arm-specific CDF:}")
  cli::cli_text(cli::col_blue(cli::style_italic("Treatment Arm")))
  print(format_dist(x$estimates$dist[1, ], x$estimates$std.error[1, ], x$estimates$ci$theta1))
  cat("\n")
  cli::cli_text(cli::col_red(cli::style_italic("Control Arm")))
  print(format_dist(x$estimates$dist[2, ], x$estimates$std.error[2, ], x$estimates$ci$theta0))
}

#' @export
print.pmf <- function(x, ...) {
  cli::cli_text("{.strong PMF Estimator}: {x$estimator}")
  cat("\n")
  cli::cli_text("{.strong Arm-specific PMF:}")
  cli::cli_text(cli::col_blue(cli::style_italic("Treatment Arm")))
  print(format_dist(x$estimates$dist[1, ], x$estimates$std.error[1, ], x$estimates$ci$theta1))
  cat("\n")
  cli::cli_text(cli::col_red(cli::style_italic("Control Arm")))
  print(format_dist(x$estimates$dist[2, ], x$estimates$std.error[2, ], x$estimates$ci$theta0))
}

#' @export
print.mannwhit <- function(x, ...) {
  cli::cli_text("{.strong Mann-Whitney Estimand:}")
  cli::cli_text(cat("     "), "{.strong Estimator}: {x$estimator}")
  cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates$theta, 2)}")
  cli::cli_text(cat("    "), "{.strong Std. error}: {round(x$estimates$std.error, 2)}")
  cli::cli_text(cat("        "), "{.strong 95% CI}: ({round(x$estimates$ci[1], 2)}, {round(x$estimates$ci[2], 2)})")
}

#' Extract RMST And Survival Probability Estimates
#'
#' @param x An object of class "rmst" or "survprob".
#'
#' @seealso \code{\link{rmst}} and \code{\link{survprob}} for creating \code{x}.
#'
#' @return A data frame containing the estimates.
#' @export
#'
#' @examples
#' \donttest{
#' surv <- survrct(Surv(time, status) ~ trt + age + sex + obstruct +
#'                    perfor + adhere + surg,
#'                 target = "trt", data = colon, coarsen = 30, estimator = "tmle")
#' est <- rmst(surv, 105:111)
#' all_estimates(est)
#' }
all_estimates <- function(x) {
  out <- do.call(
    "rbind",
    lapply(x$estimates[-which(names(x$estimates) %in%
                                c("mbcv_theta", "mbcv_treatment", "mbcv_control"))],
           function(x) {
             data.frame(treatment = x$arm1,
                        treatment.conf.low = x$arm1.conf.low,
                        treatment.conf.high = x$arm1.conf.high,
                        treatment.unif.low = x$arm1.unif.low,
                        treatment.unif.high = x$arm1.unif.high,
                        control = x$arm0,
                        control.conf.low = x$arm0.conf.low,
                        control.conf.high = x$arm0.conf.high,
                        control.unif.low = x$arm0.unif.low,
                        control.unif.high = x$arm0.unif.high,
                        theta = x$theta,
                        theta.conf.low = x$theta.conf.low,
                        theta.conf.high = x$theta.conf.high,
                        theta.unif.low = x$theta.unif.low,
                        theta.unif.high = x$theta.unif.high)
           }))
  out$horizon <- x$horizon
  out[, c(16, 1:15)]
}

format_est <- function(x) {
  x <- all_estimates(x)
  x$`treatment` <- format_digits(x$treatment, 2)
  x$`control` <- format_digits(x$control, 2)
  x$theta <- format_digits(x$theta, 2)
  x$`point-wise 95% CI` <- paste0("(", paste(format_digits(x$theta.conf.low, 2), "to", format_digits(x$theta.conf.high, 2)), ")")
  x$`uniform 95% CI` <- paste0("(", paste(format_digits(x$theta.unif.low, 2), "to", format_digits(x$theta.unif.high, 2)), ")")
  x[, c(1:2, 7, 12, 17:18)]
}

format_dist <- function(dist, std.error, ci) {
  out <- data.frame(Estimate = dist,
                    std.error = std.error,
                    ci = paste0("(", paste(format_digits(ci[, 1], 2), "to", format_digits(ci[, 2], 2)), ")"))
  names(out) <- c("Estimate", "Std. error", "95% CI")
  out
}
