metarf <- ordinalrct(Y ~ A + age + EnrollmentNIHSStotal, "A", mistie, "tmle",
                     algo = "earth", crossfit = FALSE)

log_or(metarf)
mannwhitney(metarf)

metaen <- ordinalrct(Y ~ A + age + EnrollmentNIHSStotal, "A", mistie, "tmle", algo = "lasso")

log_or(metaen)
mannwhitney(meta)

meta <- ordinalrct(Y ~ A + age + EnrollmentNIHSStotal, "A", mistie, "tmle", algo = "glm")

log_or(meta)
mannwhitney(meta)
