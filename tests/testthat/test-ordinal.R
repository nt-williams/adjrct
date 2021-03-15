meta <- ordinalrct(Y ~ A + age + EnrollmentNIHSStotal, "A", mistie, "tmle", algo = "rf")

log_or(meta)
mannwhitney(meta)

meta <- ordinalrct(Y ~ A + age + EnrollmentNIHSStotal, "A", mistie, "tmle", algo = "lasso")

log_or(meta)
mannwhitney(meta)

meta <- ordinalrct(Y ~ A + age + EnrollmentNIHSStotal, "A", mistie, "tmle", algo = "glm")

log_or(meta)
mannwhitney(meta)
