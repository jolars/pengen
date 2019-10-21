#' Fit a prototype for a penalized, generalized linear model
#'
#' @param x features
#' @param y response
#' @param family family object
#' @param penalty penalty objet
#' @param solver which solver to use
#' @param intercept whether to fit an intercept (bias)
#' @param standardize whether to standardize features
#' @param sigma regularization strength for SLOPE-type penalties
#' @param lambda regularization strength for LASSO and SLOPE
#' @param fdr target false discovery rate
#' @param n_penalties number of penalty terms to apply
#' @param lambda_min_ratio smallest size of lambda as a fraction of
#'   lambda_max
#'
#' @return A list of results.
#' @export
model <- function(x,
                  y,
                  family = "gaussian",
                  penalty = "slope",
                  sigma = NULL,
                  lambda = NULL,
                  fdr = 0.2,
                  n_penalties = 100,
                  lambda_min_ratio = ifelse(NROW(x) < NCOL(x), 0.01, 0.0001),
                  solver = "fista",
                  intercept = TRUE,
                  standardize = TRUE) {

  fit_intercept <- isTRUE(intercept)
  standardize <- isTRUE(standardize)

  family <- switch(match.arg(family),
                   gaussian = Gaussian())

  x <- as.matrix(x)
  y <- as.matrix(y)

  n <- NROW(x)
  p <- NCOL(x)
  m <- NCOL(y)

  # collect response and variable names
  response_names <- colnames(y)
  variable_names <- colnames(x)

  if (is.null(variable_names))
    variable_names <- paste0("V", seq_len(p))
  if (is.null(response_names))
    response_names <- paste0("y", seq_len(m))

  tmp <- preprocessResponse(family, y)

  y <- tmp$y
  y_scale <- tmp$y_scale
  y_center <- tmp$y_center

  tmp <- preprocessFeatures(x, standardize)

  x <- tmp$x
  x_center <- tmp$x_center
  x_scale <- tmp$x_scale

  penalty <- switch(
    match.arg(penalty),
    slope = SLOPE(x = x,
                  y = y,
                  y_scale = y_scale,
                  lambda = lambda,
                  sigma = sigma,
                  lambda_min_ratio = lambda_min_ratio,
                  n_penalties = n_penalties,
                  fdr = fdr,
                  family = family)
  )

  sigma <- penalty$sigma
  lambda <- penalty$lambda

  n_penalties <- NROW(lambda)

  intercepts <- array(0, c(1, m, n_penalties))
  betas <- array(0, c(p, m, n_penalties))

  diagnostics <- vector("list", n_penalties)

  L <- lipschitzConstant(family, x, fit_intercept)

  for (i in seq_len(n_penalties)) {
    if (i == 1) {
      beta <- matrix(0, p, m)
      intercept <- matrix(0, 1, m)
    } else {
      beta <- betas[, , i-1]
      intercept <- intercepts[, , i-1]
    }

    result <- fista(x = x,
                    y = y,
                    lambda = lambda[i, ],
                    family = family,
                    penalty = penalty,
                    intercept = intercept,
                    beta = beta,
                    standardize = standardize,
                    fit_intercept = fit_intercept,
                    L = L,
                    eta = 2,
                    tol_rel_gap = 1e-6,
                    tol_infeas = 1e-6,
                    max_passes = 1e3)

    betas[, , i] <- result$beta
    intercepts[, , i] <- result$intercept
    diagnostics[[i]] <- result$diagnostics
  }

  tmp <- unstandardize(intercepts,
                       betas,
                       x,
                       y,
                       fit_intercept,
                       x_center,
                       x_scale,
                       y_center,
                       y_scale)

  betas <- tmp$betas
  intercepts <- tmp$intercepts

  out <- list(
    coefficients = betas,
    intercept = intercepts,
    diagnostics = diagnostics,
    lambda = lambda,
    sigma = sigma
  )

  structure(out, class = "Pengen")
}
