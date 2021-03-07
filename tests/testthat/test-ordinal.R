df <- mistie

df$Y <- ordered(df$Y)
meta <- ordinalrct(Y ~ A + age, "A", df, "tmle")
log_or(meta)
cdf(meta)

drord::drord(out = as.numeric(df$Y), treat = df$A, covar = df,
             out_form = "age", treat_form = "age", param = "log_odds",
             ci = "wald", stratify = FALSE)$cdf$est

data(covid19, package = "drord")
