df <- mistie

df$Y <- ordered(df$Y)
meta <- ordinalrct(Y ~ A + age, "A", df, "unadjusted")
log_or(meta)
cdf(meta)
