infeasibility <- function(object, ...)
  UseMethod("infeasibility")

infeasibility.SLOPE <- function(object, gradient, lambda) {
  sorted_gradient <- sort(abs(gradient), decreasing = TRUE)
  max(pmax.int(cumsum(sorted_gradient - lambda)), 0)
}
