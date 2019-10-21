SLOPE <- function(x,
                  y,
                  y_scale,
                  lambda = c("gaussian", "bhq"),
                  sigma = NULL,
                  lambda_min_ratio = NULL,
                  n_penalties = 100,
                  fdr = 0.2,
                  family) {

  n <- NROW(x)
  p <- NCOL(x)

  if (is.null(lambda_min_ratio))
    lambda_min_ratio <- if (n < p) 0.01 else 0.0001

  if (is.null(lambda))
    lambda <- "gaussian"

  # noise estimate
  if (!is.null(sigma)) {
    stopifnot(length(sigma) > 0, all(sigma >= 0), is.finite(sigma))
  }

  # regularization strength
  if (is.character(lambda)) {
    lambda_type <- match.arg(lambda)

    if (lambda_type %in% c("bhq", "gaussian")) {
      q <- 1:p * fdr/(2*p)
      lambda <- stats::qnorm(1 - q)

      if (lambda_type == "gaussian" && p > 1) {
        sum_sq <- 0
        for (i in 2:p) {
          sum_sq <- sum_sq + lambda[i - 1]^2
          w <- max(1, n - i)
          lambda[i] <- lambda[i]*sqrt(1 + sum_sq/w)
        }
      }

      # ensure non-increasing lambdas
      lambda[which.min(lambda):p] <- min(lambda)

    }

    if (is.null(sigma)) {
      lambda_max <- lambdaMax(family, x, y, y_scale)

      sigma <- exp(seq(log(lambda_max/mean(lambda)),
                       log(lambda_max/mean(lambda)*lambda_min_ratio),
                       length.out = n_penalties))

      lambda <- matrix(rep(lambda, each = n_penalties, times = p), n_penalties, p)
      lambda <- lambda*sigma
    }
  } else {
    lambda <- as.matrix(lambda)
  }

  if (NCOL(lambda) != p)
    stop("lambda sequence must be as long as there are variables")

  if (any(lambda < 0))
    stop("lambda sequence cannot contain negative values")

  structure(list(name = "slope",
                 parameters = c("sigma", "fdr"),
                 sigma = sigma,
                 fdr = fdr,
                 lambda = lambda),
            class = c("SLOPE", "Penalty"))
}
