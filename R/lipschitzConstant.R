lipschitzConstant <- function(object, ...) {
  UseMethod("lipschitzConstant", object)
}

lipschitzConstant.Gaussian <- function(object,
                                       x,
                                       fit_intercept,
                                       ...) {
  norm(x^2, "I") + fit_intercept
}

lipschitzConstant.Binomial <- function(object,
                                       x,
                                       fit_intercept) {
  0.25 * (norm(x^2, "I") + fit_intercept)
}
