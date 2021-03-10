meta <- ordinalrct(Y ~ A + age, "A", mistie, "tmle", lasso = TRUE)

log_or(meta)
cdf(meta)
pmf(meta)
mannwhitney(meta)
