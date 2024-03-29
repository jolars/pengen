preprocessFeatures <- function(x, standardize) {
  p <- NCOL(x)

  if (standardize) {
    x_center <- colMeans(x)
    x_scale <- matrixStats::colSds(x)

    x_scale[x_scale == 0] <- 1

    x <- sweep(x, 2, x_center, "-")
    x <- sweep(x, 2, x_scale, "/")

    # x_center <- Matrix::colMeans(x)
    #
    # if (!inherits(x, "sparseMatrix")) {
    #   x <- sweep(x, 2, x_center)
    #   x_scale <- apply(x, 2, norm, type = "2")
    #   x_scale[x_scale == 0] <- 1
    #   x <- sweep(x, 2, x_scale, "/")
    # } else {
    #   x_scale <- standardizedSparseColNorms(x, x_center)
    #   x_scale[x_scale == 0] <- 1
    #   for (j in seq_len(p))
    #     x[, j] <- x[, j]/x_scale[j]
    # }
  } else {
    x_center <- rep.int(0, p)
    x_scale  <- rep.int(1, p)
  }

  list(x = x,
       x_center = x_center,
       x_scale = x_scale)
}
