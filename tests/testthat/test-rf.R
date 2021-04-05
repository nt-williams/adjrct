surv <- survrct(Surv(time, status) ~ trt + age + sex, target = "trt", data = colon,
                coarsen = 30, estimator = "tmle", algo = "earth")

rmst(surv, 60)
survprob(surv, 60)

surv1 <- survrct(Surv(time, status) ~ trt + age + sex, target = "trt", data = colon,
                coarsen = 30, estimator = "tmle", algo = "lasso")

rmst(surv1, 60)
survprob(surv, 60)
