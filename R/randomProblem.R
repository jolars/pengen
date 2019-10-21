#' Randomly generated problem data
#'
#' @param n numer of data points
#' @param p number of features
#' @param q proportion of nonzero features
#' @param n_groups number of groups (for group penalties)
#' @param density density of features
#' @param amplitude scale of coefficients
#' @param sigma variability in model
#' @param response type of response
#'
#' @return A list of the data along with some other stuff
#' @export
randomProblem <- function(n = 1000,
                          p = 100,
                          q = 0.2,
                          n_groups = NULL,
                          density = 1,
                          amplitude = 3,
                          sigma = 1,
                          response = c("gaussian", "binomial")) {
  if (density == 1) {
    x <- matrix(stats::rnorm(n*p), n)
  } else {
    x <- Matrix::rsparsematrix(n, p, density)
  }

  if (!is.null(n_groups)) {
    groups <- rep(seq_len(n_groups), each = ceiling(p/n_groups),
                  length.out = p)
    nonzero <- which(groups %in% seq_len(max(floor(n_groups*q), 1)))
  } else {
    groups <- NA
    nonzero <- sample(p, max(floor(q*p), 1))
  }

  beta <- amplitude * (1:p %in% nonzero)

  y <- switch(match.arg(response),
              gaussian = x %*% beta + stats::rnorm(n, sd = sigma),
              binomial = {
                y <- x %*% beta + stats::rnorm(n, sd = sigma)
                (sign(y) + 1)/2
              })

  dimnames(x) <- list(seq_len(nrow(x)),
                      paste0("V", seq_len(ncol(x))))

  list(x = x,
       y = as.double(y),
       beta = beta,
       groups = groups,
       nonzero = nonzero,
       q = q)
}
