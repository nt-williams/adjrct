x <- survrct(Surv(time, status) ~ trt + age + sex + obstruct + perfor + adhere + surg,
             target = "trt", data = colon, coarsen = 30, estimator = "tmle")

rmst(x, 20)
survprob(x, 20)

y <- survrct(Surv(time, status) ~ trt + age + sex + obstruct + perfor + adhere + surg,
             target = "trt", data = colon, coarsen = 30, estimator = "aipw")

rmst(y, 20)
survprob(y, 20)

z <- survrct(Surv(time, status) ~ trt + age + sex + obstruct + perfor + adhere + surg,
             target = "trt", data = colon, coarsen = 30, estimator = "km")

rmst(z, 20)
survprob(z, 20)
