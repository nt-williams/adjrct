
context("Various checks")

test_that("fixes NA coefficients", {
  coefs <- c(1.2, NA_real_)
  expect_equal(check_na_coef(coefs), c(1.2, 0))
})
