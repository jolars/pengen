preprocessResponse <- function(object, ...) {
  UseMethod("preprocessResponse", object)
}

preprocessResponse.Gaussian <- function(object, y, ...) {
  y <- as.numeric(y)

  if (NCOL(y) > 1)
    stop("response for Gaussian regression must be one-dimensional.")

  y_center <- mean(y)
  y_scale  <- 1

  y <- as.matrix(y - y_center)

  list(y = y,
       y_center = y_center,
       y_scale = y_scale,
       n_classes = 1L,
       class_names = NA_character_)
}
