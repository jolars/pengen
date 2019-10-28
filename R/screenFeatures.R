#' Screen features
#'
#' @param family family
#' @param penalty penalty
#' @param x features
#' @param y response
#' @param beta coefficients from last fit
#' @param intercept intercept from last fit
#' @param method type of screening
#'
#' @return A logical vector indicating whether to drop a feature or not
#' @export
screenFeatures <- function(family,
                           penalty,
                           lambda,
                           lambda_prev,
                           x,
                           y,
                           beta_prev,
                           intercept_prev,
                           method = c("none", "strong", "safe")) {
  method <- match.arg(method)

  switch(
    method,
    none = {
      rep(FALSE, NCOL(x))
    },
    strong = {
      linear_predictor <- x %*% beta_prev + intercept_prev
      grad <- gradient(family, linear_predictor, y, x)
      ord <- order(abs(grad), decreasing = TRUE)

      lh <- abs(grad)
      rh <- 2*lambda[ord] - lambda_prev[ord]

      lh < rh
    },

    safe = {
      linear_predictor <- x %*% beta_prev + intercept_prev
      pseudo_grad <- pseudoGradient(family, linear_predictor, y)
      grad <- crossprod(x, pseudo_grad)

      ord <- order(abs(grad), decreasing = TRUE)

      lh <- abs(grad)

      rh <- lambda[ord] -
        apply(x, 2, norm, "2") *
        norm(pseudo_grad, "2") *
        (lambda_prev[ord] - lambda[ord])/lambda_prev[ord]

      lh < rh
    }
  )
}
