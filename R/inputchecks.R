##' Checks to make sure the dimensions of the matrices conform
##'
##' @title Dimension check
##' @param mat_1 the U matrix
##' @param mat_2 the V matrix
##' @param dims the output of \code{dim(X)}
##' @return TRUE if it works, will error if it doesn't.
##' @keywords internal
##' @noRd
dimcheck_stop <- function(mat_1, mat_2, dims) {
  if (!symm_check(mat_1)) return("U must be symmetric.")
  if (!symm_check(mat_2)) return("V must be symmetric.")
  if (!(dims[1] == dim(mat_1)[2] && dim(mat_1)[1] == dim(mat_1)[2] &&
    dims[2] == dim(mat_2)[1] && dim(mat_2)[1] == dim(mat_2)[2])) {
      return("Non-conforming dimensions: U and V must have compatible dimensions with x.")
  }
  "ok"
}

##' Validate that the input is numeric
##'
##' Pass the arguments that are required to be numeric.
##' @title All numeric
##' @param ... The arguments that must all be numeric
##' @return TRUE if all are numeric. Otherwise, FALSE.
##' @keywords internal
##' @noRd
allnumeric_stop <- function(...) {
  args <- list(...)
  all(sapply(args, is.numeric, simplify = TRUE))
}
