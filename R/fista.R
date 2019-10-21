#' FISTA
#'
#' @param x feature matrix
#' @param y response
#' @param lambda regularization strength
#' @param family family object
#' @param penalty penalty object
#' @param intercept intercept
#' @param beta coefficients
#' @param standardize whether data is standardized
#' @param fit_intercept whether to fit intercept
#' @param L Lipschitz constant
#' @param eta step size in line search
#' @param tol_rel_gap relative tolerance for duality gap check
#' @param tol_infeas tolerance for infeasibility check
#' @param max_passes max iterations for algorithm
#'
#' @return Coefficients and diagnostics.
#' @export
fista <- function(x,
                  y,
                  lambda,
                  family,
                  penalty,
                  intercept,
                  beta,
                  standardize,
                  fit_intercept,
                  L = 1,
                  eta = 2,
                  tol_rel_gap = 1e-6,
                  tol_infeas = 1e-6,
                  max_passes = 1e3) {

  n <- NROW(x)
  p <- NCOL(x)
  m <- NCOL(y)

  intercept_tilde_old <- intercept_tilde <- intercept
  beta_tilde_old <- beta_tilde <- beta

  pseudo_gradient <- gradient <- matrix(0, p, m)

  gradient_intercept <- double(m)

  t <- 1

  accepted <- FALSE

  # diagnostics
  primals <- duals <- time <- infeasibilities <-
    rep(NA, max_passes)
  start_time <- proc.time()[3]

  linear_predictor <- x %*% beta + as.double(intercept)

  for (i in seq_len(max_passes)) {

    f <- primal(family, linear_predictor, y)
    pseudo_gradient <- pseudoGradient(family, linear_predictor, y)
    gradient <- crossprod(x, pseudo_gradient)

    if (fit_intercept)
      intercept_gradient <- mean(pseudo_gradient)

    primal_obj <- f + primal(penalty, beta, lambda)
    dual_obj <- dual(family, linear_predictor, y)
    infeasibility <- infeasibility(penalty, gradient, lambda)

    # store diagnostics
    time[i] <- proc.time()[3] - start_time[3]
    primals[i] <- primal_obj
    duals[i] <- dual_obj
    infeasibilities[i] <- infeasibility

    # check for convergence
    # TODO(JL): put this in a separate function
    infeasibility_check <- infeasibility <= tol_infeas
    duality_gap_check <-
      abs(primal_obj - dual_obj)/max(1, primal_obj) < tol_rel_gap
    accepted <- infeasibility_check && duality_gap_check

    if (accepted)
      break

    beta_tilde_old <- beta_tilde

    if (fit_intercept)
      intercept_tilde_old <- intercept_tilde

    f_old <- f
    t_old <- t

    # lipschitz search
    repeat {
      beta_tilde <- prox(penalty, beta - gradient/L, lambda, 1/L)

      d <- beta_tilde - beta

      if (fit_intercept)
        intercept_tilde <- intercept - intercept_gradient/L

      linear_predictor <- x %*% beta_tilde + as.double(intercept_tilde)

      f <- primal(family, linear_predictor, y)
      q <- f_old + crossprod(d, gradient) + 0.5 * L * sum(d^2)

      if (any(diag(q) >= f*(1 - 1e-12))) {
        break
      } else {
        L <- L*eta
      }
    }

    # FISTA step
    t <- 0.5 * (1 + sqrt(1 + 4*t_old^2))
    beta <- beta_tilde + (t_old - 1)/t * (beta_tilde - beta_tilde_old)

    if (fit_intercept) {
      intercept <- intercept_tilde +
        (t_old - 1)/t * (intercept_tilde - intercept_tilde_old)
    }

    linear_predictor <- x %*% beta + as.double(intercept)
  }


  diagnostics <- list(time = time[seq_len(i)],
                      primal = duals[seq_len(i)],
                      dual = primals[seq_len(i)],
                      infeasibility = infeasibilities[seq_len(i)])

  list(beta = beta,
       intercept = intercept,
       passes = i,
       diagnostics = diagnostics)
}
