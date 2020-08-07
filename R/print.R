
#' @export
print.rmst <- function(x, ...) {
  cli::cli_text("{.strong RMST Estimator}: {x$estimator}")
  if (length(x$horizon) > 1) {
    cat("\n")
    cat("           Confidence level: 95%\n")
    cat("     Multiplier Bootstrap C:", "\n")
    cat(" Test of no effect, p-value:\n")
    cat("      First 6 time horizons:\n")
    cat("\n")
    print(head(all_rmst(x)))
    cli::cli_text(cli::col_red("Access all estimates with `all_rmst()`"))
  } else {
    cli::cli_text(cat("  "), "{.strong Time horizon}: {x$horizon}")
    cat("\n")
    cli::cli_text("{.strong Arm-specific RMST:}")
    cli::cli_text(cli::col_blue(cli::style_italic("Treatment Arm")))
    cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates[[1]]$rmst1, 2)}")
    cli::cli_text(cat("    "), "{.strong Std. error}: ")
    cli::cli_text(cat("        "), "{.strong 95% CI}: ")
    cli::cli_text(cli::col_red(cli::style_italic("Control Arm")))
    cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates[[1]]$rmst0, 2)}")
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
all_rmst <- function(x) {
  out <- do.call("rbind", lapply(x$estimates, function(x)
    data.frame(
      treatment       = x$rmst1,
      control         = x$rmst0,
      theta           = x$theta,
      theta.conf.low  = x$theta.conf.low,
      theta.conf.high = x$theta.conf.high
    )))
  out$horizon <- 1 + (1:length(x$estimates))
  out[, c(6, 1:5)]
}
