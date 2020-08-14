
check_name_conflicts <- function(data) {


}

check_na_coef <- function(coefs) {
  coefs[is.na(coefs)] <- 0
  return(coefs)
}

check_sl3 <- function() {
  if (length(find.package("sl3", quiet = TRUE)) == 0) {
    FALSE
  } else {
    TRUE
  }
}
