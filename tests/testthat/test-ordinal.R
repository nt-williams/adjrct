library(devtools)

load_all()

df <- mistie

df$Y <- ordered(df$Y)
meta <- ordinalrct(Y ~ A + age + EnrollmentNIHSStotal, "A", df, "tmle", lasso = TRUE)

log_or(meta)
cdf(meta)
pmf(meta)
mannwhitney(meta)

x <- cdf(meta)$estimates




cli::cli_text("{.strong Arm-specific CDF:}")
cli::cli_text(cli::col_blue(cli::style_italic("Treatment Arm")))
print(format_dist(x$dist[1, ], x$std.error[1, ], x$ci$theta1))
cat("\n")
cli::cli_text(cli::col_red(cli::style_italic("Control Arm")))
print(format_dist(x$dist[2, ], x$std.error[2, ], x$ci$theta0))


print.cdf(cdf(meta))
