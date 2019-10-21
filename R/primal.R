primal <- function(object, ...) {
  UseMethod("primal", object)
}

primal.Gaussian <- function(object, linear_predictor, y, ...) {
  0.5 * norm(linear_predictor - y, "2")^2
}

primal.SLOPE <- function(object, beta, lambda, ...) {
  sum(lambda * sort(abs(beta), decreasing = TRUE))
}
