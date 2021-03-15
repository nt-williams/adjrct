surv <- survrct(Surv(time, status) ~ trt + age + sex, target = "trt", data = colon,
                coarsen = 30, estimator = "tmle", algo = "rf")

rmst(surv, 60)
survprob(surv, 60)
