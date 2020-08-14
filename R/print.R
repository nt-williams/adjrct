
#' @export
print.rmst <- function(x, ...) {
  cli::cli_text("{.strong RMST Estimator}: {x$estimator}")
  if (length(x$horizon) > 1) {
    cat("\n")
    cat("           Confidence level: 95%\n")
    cat("     Multiplier Bootstrap C:", x$estimates$mbcv, "\n")
    cat(" Test of no effect, p-value:\n")
    cat("          First 6 estimates:\n")
    cat("\n")
    print(head(pretty_print(x)))
    cli::cli_text(cli::col_red("Access all estimates with `all_estimates()`"))
  } else {
    cli::cli_text(cat("  "), "{.strong Time horizon}: {x$horizon}")
    cat("\n")
    cli::cli_text("{.strong Arm-specific RMST:}")
    cli::cli_text(cli::col_blue(cli::style_italic("Treatment Arm")))
    cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates[[1]]$arm1, 2)}")
    cli::cli_text(cat("    "), "{.strong Std. error}: ")
    cli::cli_text(cat("        "), "{.strong 95% CI}: ")
    cli::cli_text(cli::col_red(cli::style_italic("Control Arm")))
    cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates[[1]]$arm0, 2)}")
    cli::cli_text(cat("    "), "{.strong Std. error}: ")
    cli::cli_text(cat("        "), "{.strong 95% CI}: ")
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
    cat("     Multiplier Bootstrap C:", x$estimates$mbcv, "\n")
    cat(" Test of no effect, p-value:\n")
    cat("          First 6 estimates:\n")
    cat("\n")
    print(head(pretty_print(x)))
    cli::cli_text(cli::col_red("Access all estimates with `all_estimates()`"))
  } else {
    cli::cli_text(cat("  "), "{.strong Time horizon}: {x$horizon}")
    cat("\n")
    cli::cli_text("{.strong Arm-specific Survival Probability:}")
    cli::cli_text(cli::col_blue(cli::style_italic("Treatment Arm")))
    cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates[[1]]$arm1, 2)}")
    cli::cli_text(cat("    "), "{.strong Std. error}: ")
    cli::cli_text(cat("        "), "{.strong 95% CI}: ")
    cli::cli_text(cli::col_red(cli::style_italic("Control Arm")))
    cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates[[1]]$arm0, 2)}")
    cli::cli_text(cat("    "), "{.strong Std. error}: ")
    cli::cli_text(cat("        "), "{.strong 95% CI}: ")
    cat("\n")
    cli::cli_text("{.strong Treatment Effect:}")
    cli::cli_text(cli::col_green(cli::style_italic("Additive effect")))
    cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates[[1]]$theta, 2)}")
    cli::cli_text(cat("    "), "{.strong Std. error}: {round(x$estimates[[1]]$std.error, 2)}")
    cli::cli_text(cat("        "), "{.strong 95% CI}: ({round(x$estimates[[1]]$theta.conf.low, 2)}, {round(x$estimates[[1]]$theta.conf.high, 2)})")
  }
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
all_estimates <- function(x) {
  out <- do.call("rbind", lapply(x$estimates[-which(names(x$estimates) == "mbcv")], function(x)
    data.frame(
      treatment       = x$arm1,
      control         = x$arm0,
      theta           = x$theta,
      theta.conf.low  = x$theta.conf.low,
      theta.conf.high = x$theta.conf.high,
      unif.conf.low   = x$unif.conf.low,
      unif.conf.high  = x$unif.conf.high
    )))
  out$horizon <- x$horizon
  out[, c(8, 1:7)]
}

pretty_print <- function(x) {
  x <- all_estimates(x)
  x$`Treatment` <- format_digits(x$treatment, 2)
  x$`Control` <- format_digits(x$control, 2)
  x$Theta <- format_digits(x$theta, 2)
  x$`Point-wise 95% CI` <- paste0("(", paste(format_digits(x$theta.conf.low, 2), "to", format_digits(x$theta.conf.high, 2)), ")")
  x$`Uniform 95% CI` <- paste0("(", paste(format_digits(x$unif.conf.low, 2), "to", format_digits(x$unif.conf.high, 2)), ")")
  x[, c(1, 9:13)]
}

