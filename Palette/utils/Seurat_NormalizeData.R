library(future)
# Generate chunk points
#
# @param dsize How big is the data being chunked
# @param csize How big should each chunk be
#
# @return A matrix where each column is a chunk, row 1 is start points, row 2 is end points
#
ChunkPoints1 <- function(dsize, csize) {
  return(vapply(
    X = 1L:ceiling(x = dsize / csize),
    FUN = function(i) {
      return(c(
        start = (csize * (i - 1L)) + 1L,
        end = min(csize * i, dsize)
      ))
    },
    FUN.VALUE = numeric(length = 2L)
  ))
}

NormalizeData1 <- function(
    object,
    normalization.method = "LogNormalize",
    scale.factor = 1e4,
    margin = 1,
    block.size = NULL,
    verbose = TRUE,
    ...
) {
  # CheckDots(...)
  if (is.null(x = normalization.method)) {
    return(object)
  }
  normalized.data <- if (nbrOfWorkers() > 1) {
    norm.function <- switch(
      EXPR = normalization.method,
      'LogNormalize' = LogNormalize1,
      'CLR' = CustomNormalize1,
      stop("Unknown normalization method: ", normalization.method)
    )
    if (normalization.method != 'CLR') {
      margin <- 2
    }
    dsize <- switch(
      EXPR = margin,
      '1' = nrow(x = object),
      '2' = ncol(x = object),
      stop("'margin' must be 1 or 2")
    )
    chunk.points <- ChunkPoints1(
      dsize = dsize,
      csize = block.size %||% ceiling(x = dsize / nbrOfWorkers())
    )
    normalized.data <- future_lapply(
      X = 1:ncol(x = chunk.points),
      FUN = function(i) {
        block <- chunk.points[, i]
        data <- if (margin == 1) {
          object[block[1]:block[2], , drop = FALSE]
        } else {
          object[, block[1]:block[2], drop = FALSE]
        }
        clr_function <- function(x) {
          return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))))
        }
        args <- list(
          data = data,
          scale.factor = scale.factor,
          verbose = FALSE,
          custom_function = clr_function, margin = margin
        )
        args <- args[names(x = formals(fun = norm.function))]
        return(do.call(
          what = norm.function,
          args = args
        ))
      }
    )
    do.call(
      what = switch(
        EXPR = margin,
        '1' = 'rbind',
        '2' = 'cbind',
        stop("'margin' must be 1 or 2")
      ),
      args = normalized.data
    )
  } else {
    switch(
      EXPR = normalization.method,
      'LogNormalize' = LogNormalize1(
        data = object,
        scale.factor = scale.factor,
        verbose = verbose
      ),
      'CLR' = CustomNormalize1(
        data = object,
        custom_function = function(x) {
          return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))))
        },
        margin = margin,
        verbose = verbose
        # across = across
      ),
      stop("Unknown normalization method: ", normalization.method)
    )
  }
  return(normalized.data)
}

LogNorm <- function(data, scale_factor, display_progress = TRUE) {
  .Call('_Seurat_LogNorm', PACKAGE = 'Seurat', data, scale_factor, display_progress)
}
#' Normalize raw data
#'
#' Normalize count data per cell and transform to log scale
#'
#' @param data Matrix with the raw count data
#' @param scale.factor Scale the data. Default is 1e4
#' @param verbose Print progress
#'
#' @return Returns a matrix with the normalize and log transformed data
#'
#' @importFrom methods as
#'
#' @export
#' @concept preprocessing
#'
#' @examples
#' mat <- matrix(data = rbinom(n = 25, size = 5, prob = 0.2), nrow = 5)
#' mat
#' mat_norm <- LogNormalize(data = mat)
#' mat_norm
#'
LogNormalize1 <- function(data, scale.factor = 1e4, verbose = TRUE) {
  if (is.data.frame(x = data)) {
    data <- as.matrix(x = data)
  }
  if (!inherits(x = data, what = 'dgCMatrix')) {
    data <- as.sparse(x = data)
  }
  # call Rcpp function to normalize
  if (verbose) {
    cat("Performing log-normalization\n", file = stderr())
  }
  norm.data <- LogNorm(data, scale_factor = scale.factor, display_progress = verbose)
  colnames(x = norm.data) <- colnames(x = data)
  rownames(x = norm.data) <- rownames(x = data)
  return(norm.data)
}



# Normalize a given data matrix
#
# Normalize a given matrix with a custom function. Essentially just a wrapper
# around apply. Used primarily in the context of CLR normalization.
#
# @param data Matrix with the raw count data
# @param custom_function A custom normalization function
# @param margin Which way to we normalize. Set 1 for rows (features) or 2 for columns (genes)
# @parm across Which way to we normalize? Choose form 'cells' or 'features'
# @param verbose Show progress bar
#
# @return Returns a matrix with the custom normalization
#
#' @importFrom Matrix t
#' @importFrom methods as
#' @importFrom pbapply pbapply
#

CustomNormalize1 <- function(data, custom_function, margin, verbose = TRUE) {
  if (is.data.frame(x = data)) {
    data <- as.matrix(x = data)
  }
  if (!inherits(x = data, what = 'dgCMatrix')) {
    data <- as.sparse(x = data)
  }
  myapply <- ifelse(test = verbose, yes = pbapply, no = apply)
  if (verbose) {
    message("Normalizing across ", c('features', 'cells')[margin])
  }
  norm.data <- myapply(
    X = data,
    MARGIN = margin,
    FUN = custom_function)
  if (margin == 1) {
    norm.data = Matrix::t(x = norm.data)
  }
  colnames(x = norm.data) <- colnames(x = data)
  rownames(x = norm.data) <- rownames(x = data)
  return(norm.data)
}


# example
#d1 <- readRDS("/Users/applexiufu/Library/CloudStorage/OneDrive-stu.hit.edu.cn/data/dataset10/dataset.rds")
#a <- NormalizeData1(d1[[1]])
#a <- LogNormalize1(d1[[1]],1e4,T)
