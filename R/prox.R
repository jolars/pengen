#' Evaluate proximal operator
#'
#' @param object penalty object, see [Penalty].
#' @param ... ignored
#'
#' @return Coefficients after evaluating prox
#' @export
prox <- function(object, ...) {
  UseMethod("prox")
}

prox.SLOPE <- function(object, coefs, lambda, shrinkage, ...) {
  slopeProx(coefs, lambda, shrinkage)
}
