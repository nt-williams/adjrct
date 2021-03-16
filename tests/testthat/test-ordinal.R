metarf <- ordinalrct(Y ~ A + age + EnrollmentNIHSStotal, "A", mistie, "tmle", algo = "rf")

log_or(metarf)
mannwhitney(meta)

metaen <- ordinalrct(Y ~ A + age + EnrollmentNIHSStotal, "A", mistie, "tmle", algo = "lasso")

log_or(metaen)
mannwhitney(meta)

meta <- ordinalrct(Y ~ A + age + EnrollmentNIHSStotal, "A", mistie, "tmle", algo = "glm")

log_or(meta)
mannwhitney(meta)
