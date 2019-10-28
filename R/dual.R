 dual <- function(object, ...) {
  UseMethod("dual", object)
}

dual.Gaussian <- function(object, linear_predictor, y, ...) {
  # r <- linear_predictor - y
  # -norm(r, "2")^2 - sum(r*y)
  0.5*sum(y^2) - 0.5*norm(y - linear_predictor, "2")^2
}
