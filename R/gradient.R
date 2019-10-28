#' Pseudo-Gradient
#'
#' Compute the gradient (before taking the inner product with the
#' feature matrix).
#'
#' @param linear_predictor linear predictor
#' @param y response
#' @param family a family, such as [Gaussian()]
#'
#' @return Pseudo gradient (needs to be multiplied with X^T)
#' @export
pseudoGradient <- function(family, linear_predictor, y) {
  UseMethod("pseudoGradient", family)
}

pseudoGradient.Gaussian <- function(family, linear_predictor, y) {
   -(y - linear_predictor)
}

#' Gradient
#'
#' @param family family
#' @param linear_predictor linear predictor
#' @param y response
#' @param x features
#'
#' @return Gradient vector
#' @export
gradient <- function(family, linear_predictor, y, x) {
  UseMethod("gradient", family)
}

gradient.Gaussian <- function(family, linear_predictor, y, x) {
  crossprod(x, pseudoGradient(family, linear_predictor, y))
}
