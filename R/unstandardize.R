#' Return coefficients to original scale
#'
#' @param intercepts intercepts
#' @param betas coefficients
#' @param x features
#' @param y response
#' @param fit_intercept whether intercept have been fitted
#' @param x_center centers of features
#' @param x_scale scale of features
#' @param y_center centers of responses
#' @param y_scale scale of responses
#'
#' @return Coefficients and intercept on original scale.
#' @export
unstandardize <- function(intercepts,
                          betas,
                          x,
                          y,
                          fit_intercept,
                          x_center,
                          x_scale,
                          y_center,
                          y_scale) {

  p <- NROW(betas)
  m <- NCOL(betas)
  n_penalties <- dim(betas)[3]

  for (k in seq_len(m)) {
    x_bar_beta_sum <- double(n_penalties)

    for (j in seq_len(p)) {
      betas[j, k, ] <- betas[j, k, ] * y_scale[k]/x_scale[j]
      x_bar_beta_sum <- x_bar_beta_sum + x_center[j] * betas[j, k, ]
    }

    if (fit_intercept)
      intercepts[, k, ] <-
        intercepts[, k, ]*y_scale[k] + y_center[k] - x_bar_beta_sum
  }

  list(intercepts = intercepts,
       betas = betas)
}
