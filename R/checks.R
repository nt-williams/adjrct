
check_na_coef <- function(coefs) {
  coefs[is.na(coefs)] <- 0
  return(coefs)
}
