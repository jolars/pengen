#' Compute the "pseudo-gradient"
#'
#' Compute the gradient (before taking the inner product with the
#' feature matrix).
#'
#' @param family a family, such as [Gaussian()]
#' @param ... ignored
#'
#' @return Return the pseudo-gradient.
#' @export
pseudoGradient <- function(family, ...) {
  UseMethod("pseudoGradient", family)
}

pseudoGradient.Gaussian <- function(family, linear_predictor, y, ...) {
   linear_predictor - y
}
