
#' @export
print.rct_rmst <- function(x, ...) {
  cli::cli_text("{.strong RMST Estimator}: {x$estimator}")
  cli::cli_text(cat("  "), "{.strong Time horizon}: {x$horizon}")
  cat("\n")
  cli::cli_text("{.strong Arm-specific RMST:}")
  cli::cli_text(cli::col_blue(cli::style_italic("Treatment Arm")))
  cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates$rmst1, 2)}")
  cli::cli_text(cat("    "), "{.strong Std. error}: ")
  cli::cli_text(cat("        "), "{.strong 95% CI}: ")
  cli::cli_text(cli::col_red(cli::style_italic("Control Arm")))
  cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates$rmst0, 2)}")
  cli::cli_text(cat("    "), "{.strong Std. error}: ")
  cli::cli_text(cat("        "), "{.strong 95% CI}: ")
  cat("\n")
  cli::cli_text("{.strong Treatment Effect:}")
  cli::cli_text(cli::col_green(cli::style_italic("Additive effect")))
  cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates$theta, 2)}")
  cli::cli_text(cat("    "), "{.strong Std. error}: {round(x$estimates$std.error, 2)}")
  cli::cli_text(cat("        "), "{.strong 95% CI}: ({round(x$estimates$theta.conf.low, 2)}, {round(x$estimates$theta.conf.high, 2)})")
}
