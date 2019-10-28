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
                  screening = c("none", "strong", "safe"),
                  solver = "fista",
                  intercept = TRUE,
                  standardize = TRUE) {

  fit_intercept <- isTRUE(intercept)
  standardize <- isTRUE(standardize)

  screening <- match.arg(screening)

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
  lambdas <- penalty$lambda

  n_penalties <- NROW(lambdas)

  intercepts <- matrix(0, n_penalties, m)
  betas <- array(0, c(p, m, n_penalties))

  diagnostics <- vector("list", n_penalties)

  L <- lipschitzConstant(family, x, fit_intercept)

  active_set <- seq_len(p)
  ever_active <- matrix(FALSE, p, n_penalties)

  for (i in seq_len(n_penalties)) {
    if (i == 1) {
      beta <- matrix(0, p, m)
      intercept <- rep(0, m)

      lambda_prev <- lambdas[i, ]
    } else {
      beta <- betas[, , i-1]
      intercept <- intercepts[i-1, ]

      lambda_prev <- lambdas[i-1, ]
    }

    lambda <- lambdas[i, ]


    # screen features to see if they can be excluded prior to fitting
    inactive <- screenFeatures(family,
                               penalty,
                               lambda,
                               lambda_prev,
                               x,
                               y,
                               beta_prev = beta,
                               intercept_prev = as.double(intercept),
                               method = screening)
    active <- !inactive

    # keep a set of predictors that has ever been active
    ever_active[, i] <- ever_active[, max(i-1, 1)] | active

    active_set <- which(active)
    ever_active_set <- which(ever_active[, i])

    if (any(ever_active)) {
      result <- fista(x = x[, ever_active_set],
                      y = y,
                      lambda = lambdas[i, seq_along(ever_active_set)],
                      family = family,
                      penalty = penalty,
                      intercept = intercept,
                      beta = as.matrix(beta[ever_active_set]),
                      standardize = standardize,
                      fit_intercept = fit_intercept,
                      L = L,
                      eta = 2,
                      tol_rel_gap = 1e-6,
                      tol_infeas = 1e-6,
                      max_passes = 1e3)

      betas[ever_active_set, , i] <- result$beta
      intercepts[i, ] <- result$intercept
      diagnostics[[i]] <- result$diagnostics
    } else {
      diagnostics[[i]] <- NULL
    }
  }

  # return coefficients on original scale
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
    lambda = lambdas,
    sigma = sigma,
    ever_active = ever_active
  )

  structure(out, class = "Pengen")
}
