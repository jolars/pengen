#' Smallest L1-regularization penalty where solution is completely sparse
#'
#' @param object an object of class `"Family"`
#' @param x feature matrix
#' @param y response
#' @param y_scale scale of response
#'
#' @return Lambda_max
#'
#' @keywords internal
lambdaMax <- function(object, x, y, y_scale) {
  UseMethod("lambdaMax")
}

#' @rdname lambdaMax
lambdaMax.Gaussian <- function(object, x, y, y_scale) {
  y_scale * max(abs(crossprod(x, y)))
}

#' @rdname lambdaMax
lambdaMax.Binomial <- function(object, x, y, y_scale) {
  # convert y from {-1, 1} to {0, 1}
  y <- (y + 1) / 2

  # standardize
  y_center <- mean(y)
  y_scale <- stats::sd(y)

  y <- (y - y_center) / y_scale

  y_scale * max(abs(crossprod(x, y)))
}
