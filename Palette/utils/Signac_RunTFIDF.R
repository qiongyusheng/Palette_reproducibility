RunTFIDF1 <- function(
    object,
    scale.factor = 1e4,
    idf = NULL,
    verbose = TRUE,
    ...
) {
  if (inherits(x = object, what = "data.frame")) {
    object <- as.matrix(x = object)
  }
  if (!inherits(x = object, what = "CsparseMatrix")) {
    object <- as(object = object, Class = "CsparseMatrix")
  }
  if (verbose) {
    message("Performing TF-IDF normalization")
  }
  npeaks <- colSums(x = object)
  
  if (any(npeaks == 0)) {
    warning("Some cells contain 0 total counts")
  }
  #tf <- object %*% Matrix::Diagonal(x = 1 / npeaks)
  tf <- tcrossprod(x = object, y = Matrix::Diagonal(x = 1 / npeaks))

  if (!is.null(x = idf)) {
    precomputed_idf <- TRUE
    if (!inherits(x = idf, what = "numeric")) {
      stop("idf parameter must be a numeric vector")
    }
    if (length(x = idf) != nrow(x = object)) {
      stop("Length of supplied IDF vector does not match",
           " number of rows in input matrix")
    }
    if (any(idf == 0)) {
      stop("Supplied IDF values cannot be zero")
    }
    if (verbose) {
      message("Using precomputed IDF vector")
    }
  } else {
    precomputed_idf <- FALSE
    rsums <- rowSums(x = object)
    if (any(rsums == 0)) {
      warning("Some features contain 0 total counts")
    }
    idf <- ncol(x = object) / rsums
  }
  
  norm.data <- Matrix::Diagonal(n = length(x = idf), x = idf) %*% tf
  slot(object = norm.data, name = "x") <- log1p(
    x = slot(object = norm.data, name = "x") * scale.factor)
  colnames(x = norm.data) <- colnames(x = object)
  rownames(x = norm.data) <- rownames(x = object)
  # set NA values to 0
  vals <- slot(object = norm.data, name = "x")
  vals[is.na(x = vals)] <- 0
  slot(object = norm.data, name = "x") <- vals
  return(norm.data)
}

