df <- mistie

df$Y <- ordered(df$Y)
meta <- ordinalrct(Y ~ A + age, "A", df, "tmle")
log_or(meta)

cdf(meta)
pmf(meta)
