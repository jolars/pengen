#' Obtain Coefficients from Fit
#'
#' This function is equivalent to simply calling `drop(object$coefficients)`.
#'
#' @param object an object of class `'Pengen'`.
#' @param ... arguments that are passed on to [stats::update()].
#' @inheritParams predict.Golem
#'
#' @return Coefficients from the model after having dropped extraneous
#'   dimensions by calling drop.
#'
#' @export
#' @examples
#' fit <- golem(mtcars$mpg, mtcars$vs, n_sigma = 1)
#' coef(fit)
coef.Pengen <- function(object,
                       simplify = TRUE,
                       ...) {
  beta <- object$coefficients

  if (simplify)
    beta <- drop(beta)

  beta
}
