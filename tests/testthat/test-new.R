x <- survrct(Surv(time, status) ~ trt + age + sex + obstruct + perfor + adhere + surg,
             target = "trt", data = colon, coarsen = 30, estimator = "tmle")

rmst(x, 20)
survprob(x, 20)
