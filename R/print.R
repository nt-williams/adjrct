#' @export
print.rmst <- function(x, ...) {
  cli::cli_text("{.strong RMST Estimator}: {x$estimator}")
  if (length(x$horizon) > 1) {
    cat("\n")
    cat("           Confidence level: 95%\n")
    cat("          Mult. Bootstrap C:", x$estimates$mbcv_theta, "\n")
    cat("\n")
    print(format_est(x))
  } else {
    cat("\n")
    cli::cli_text("{.strong Marginal RMST:} E(min[T, {x$horizon}] | A = a)")
    cli::cli_text(cli::col_blue(cli::style_italic("Treatment Arm")))
    cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates[[1]]$arm1, 2)}")
    cli::cli_text(cat("    "), "{.strong Std. error}: {round(x$estimates[[1]]$arm1.std.error, 2)}")
    cli::cli_text(cat("        "), "{.strong 95% CI}: ({round(x$estimates[[1]]$arm1.conf.low, 2)}, {round(x$estimates[[1]]$arm1.conf.high, 2)})")
    cli::cli_text(cli::col_red(cli::style_italic("Control Arm")))
    cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates[[1]]$arm0, 2)}")
    cli::cli_text(cat("    "), "{.strong Std. error}: {round(x$estimates[[1]]$arm0.std.error, 2)}")
    cli::cli_text(cat("        "), "{.strong 95% CI}: ({round(x$estimates[[1]]$arm0.conf.low, 2)}, {round(x$estimates[[1]]$arm0.conf.high, 2)})")
    cat("\n")
    cli::cli_text("{.strong Treatment Effect:} E(min[T, {x$horizon}] | A = 1) - E(min[T, {x$horizon}] | A = 0)")
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
    print(format_est(x))
  } else {
    cat("\n")
    cli::cli_text("{.strong Marginal Survival Probability:} Pr(T > {x$horizon} | A = a)")
    cli::cli_text(cli::col_blue(cli::style_italic("Treatment Arm")))
    cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates[[1]]$arm1, 2)}")
    cli::cli_text(cat("    "), "{.strong Std. error}: {round(x$estimates[[1]]$arm1.std.error, 2)}")
    cli::cli_text(cat("        "), "{.strong 95% CI}: ({round(x$estimates[[1]]$arm1.conf.low, 2)}, {round(x$estimates[[1]]$arm1.conf.high, 2)})")
    cli::cli_text(cli::col_red(cli::style_italic("Control Arm")))
    cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates[[1]]$arm0, 2)}")
    cli::cli_text(cat("    "), "{.strong Std. error}: {round(x$estimates[[1]]$arm0.std.error, 2)}")
    cli::cli_text(cat("        "), "{.strong 95% CI}: ({round(x$estimates[[1]]$arm0.conf.low, 2)}, {round(x$estimates[[1]]$arm0.conf.high, 2)})")
    cat("\n")
    cli::cli_text("{.strong Treatment Effect:} Pr(T > {x$horizon} | A = 1) - Pr(T > {x$horizon} | A = 0)")
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
    cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates$lo1, 2)}")
    cli::cli_text(cat("    "), "{.strong Std. error}: {round(x$estimates$lo1.std.error, 2)}")
    cli::cli_text(cat("        "), "{.strong 95% CI}: ({round(x$estimates$lo1.ci[1], 2)}, {round(x$estimates$lo1.ci[2], 2)})")
    cli::cli_text(cli::col_red(cli::style_italic("Control Arm")))
    cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates$lo0, 2)}")
    cli::cli_text(cat("    "), "{.strong Std. error}: {round(x$estimates$lo0.std.error, 2)}")
    cli::cli_text(cat("        "), "{.strong 95% CI}: ({round(x$estimates$lo0.ci[1], 2)}, {round(x$estimates$lo0.ci[2], 2)})")
    cat("\n")
    cli::cli_text("{.strong Average log odds ratio:}")
    cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates$lor, 2)}")
    cli::cli_text(cat("    "), "{.strong Std. error}: {round(x$estimates$std.error, 2)}")
    cli::cli_text(cat("        "), "{.strong 95% CI}: ({round(x$estimates$ci[1], 2)}, {round(x$estimates$ci[2], 2)})")
}

#' @export
print.cdf <- function(x, ...) {
  cli::cli_text("{.strong CDF Estimator}: {x$estimator}")
  cat("\n")
  cli::cli_text("{.strong Arm-specific CDF:} Pr(K <= k | A = a)")
  cli::cli_text(cli::col_blue(cli::style_italic("Treatment Arm")))
  print(format_dist(x$estimates$dist[1, ], x$estimates$std.error[1, ], x$estimates$ci$theta1, x$estimates$ci$unif1, levels = x$levels, "cdf"))
  cat("\n")
  cli::cli_text(cli::col_red(cli::style_italic("Control Arm")))
  print(format_dist(x$estimates$dist[2, ], x$estimates$std.error[2, ], x$estimates$ci$theta0, x$estimates$ci$unif0, levels = x$levels, "cdf"))
}

#' @export
print.pmf <- function(x, ...) {
  cli::cli_text("{.strong PMF Estimator}: {x$estimator}")
  cat("\n")
  cli::cli_text("{.strong Arm-specific PMF:} Pr(K = k | A = a)")
  cli::cli_text(cli::col_blue(cli::style_italic("Treatment Arm")))
  print(format_dist(x$estimates$dist[1, ], x$estimates$std.error[1, ], x$estimates$ci$theta1, x$estimates$ci$unif1, levels = x$levels, "pmf"))
  cat("\n")
  cli::cli_text(cli::col_red(cli::style_italic("Control Arm")))
  print(format_dist(x$estimates$dist[2, ], x$estimates$std.error[2, ], x$estimates$ci$theta0, x$estimates$ci$unif0, levels = x$levels, "pmf"))
}

#' @export
print.mannwhit <- function(x, ...) {
  cli::cli_text("{.strong Mann-Whitney Estimand}")
  cat("\n")
  cli::cli_text(cat("     "), "{.strong Estimator}: {x$estimator}")
  cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$estimates$mann.whitney, 2)}")
  cli::cli_text(cat("    "), "{.strong Std. error}: {round(x$estimates$std.error, 2)}")
  cli::cli_text(cat("        "), "{.strong 95% CI}: ({round(x$estimates$ci[1], 2)}, {round(x$estimates$ci[2], 2)})")
  cat("\n")
  cat(cli::cli_text(cli::col_green("The probability that a randomly drawn patient from the treatment arm has a greater outcome level than a randomly drawn patient from the control arm, with ties broken at random.")))
}

format_est <- function(x) {
  tidied <- tidy(x)
  if (class(x) == "rmst") {
    treatment <- format_digits(tidied$trt.rmst, 2)
    control <- format_digits(tidied$control.rmst, 2)
  } else {
    treatment <- format_digits(tidied$trt.survprob, 2)
    control <- format_digits(tidied$control.survprob, 2)
  }
  data.frame(
    horizon = x$horizon,
    treatment = treatment,
    control = control,
    theta = format_digits(tidied$theta, 2),
    `point-wise 95% CI` =
      paste0("(", paste(
        format_digits(tidied$theta.conf.low, 2),
        "to",
        format_digits(tidied$theta.conf.high, 2)
      ), ")"),
    `uniform 95% CI` = paste0("(", paste(
      format_digits(tidied$theta.unif.low, 2),
      "to",
      format_digits(tidied$theta.unif.high, 2)
    ), ")"),
    check.names = FALSE
  )
}

format_dist <- function(dist, std.error, ci, unif, levels, type = c("cdf", "pmf")) {
  if (match.arg(type) == "cdf") {
    out <- data.frame(k = levels,
                      Estimate = c(format_digits(dist, 3), "1.000"),
                      std.error = c(format_digits(std.error, 3), "-"),
                      ci = c(paste0("(", paste(format_digits(ci[, 1], 2), "to", format_digits(ci[, 2], 2)), ")"), "-"),
                      unif = c(paste0("(", paste(format_digits(unif[, 1], 2), "to", format_digits(unif[, 2], 2)), ")"), "-"))
  } else {
    out <- data.frame(k = levels,
                      Estimate = format_digits(dist, 3),
                      std.error = format_digits(std.error, 3),
                      ci = paste0("(", paste(format_digits(ci[, 1], 2), "to", format_digits(ci[, 2], 2)), ")"),
                      unif = paste0("(", paste(format_digits(unif[, 1], 2), "to", format_digits(unif[, 2], 2)), ")"))
  }
  names(out) <- c("k", "Estimate", "Std. error", "95% CI", "Uniform 95% CI")
  out
}
