SparseRowVar2 <- function(mat, mu, display_progress) {
  .Call('_Seurat_SparseRowVar2', PACKAGE = 'Seurat', mat, mu, display_progress)
}
SparseRowVarStd <- function(mat, mu, sd, vmax, display_progress) {
  .Call('_Seurat_SparseRowVarStd', PACKAGE = 'Seurat', mat, mu, sd, vmax, display_progress)
}

FindVariableFeatures1 <- function(
    object,
    loess.span = 0.3,
    clip.max = 'auto',
    nfeatures = 2000,
    verbose = TRUE,
    ...
) {
  
  if (!inherits(x = object, 'Matrix')) {
    object <- as(object = as.matrix(x = object), Class = 'Matrix')
  }
  if (!inherits(x = object, what = 'dgCMatrix')) {
    object <- as.sparse(x = object)
  }
  
  if (clip.max == 'auto') {
    clip.max <- sqrt(x = ncol(x = object))
  }
  hvf.info <- data.frame(mean = rowMeans(x = object))
  hvf.info$variance <- SparseRowVar2(
    mat = object,
    mu = hvf.info$mean,
    display_progress = verbose
  )
  hvf.info$variance.expected <- 0
  hvf.info$variance.standardized <- 0
  not.const <- hvf.info$variance > 0
  fit <- loess(
    formula = log10(x = variance) ~ log10(x = mean),
    data = hvf.info[not.const, ],
    span = loess.span
  )
  hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted
  # use c function to get variance after feature standardization
  hvf.info$variance.standardized <- SparseRowVarStd(
    mat = object,
    mu = hvf.info$mean,
    sd = sqrt(hvf.info$variance.expected),
    vmax = clip.max,
    display_progress = verbose
  )
  colnames(x = hvf.info) <- paste0('vst.', colnames(x = hvf.info))
  
  
  hvf.info <- hvf.info[which(x = hvf.info[, 1, drop = TRUE] != 0), ]
  hvf.info <- hvf.info[order(hvf.info$vst.variance.standardized, decreasing = TRUE), , drop = FALSE]
  
  top.features <- head(x = rownames(x = hvf.info), n = nfeatures)
  
  return(top.features)
}


# example
# d1 <- readRDS("/Users/applexiufu/Library/CloudStorage/OneDrive-stu.hit.edu.cn/data/dataset10/dataset.rds")
# a <- FindVariableFeatures1(d1[[1]])


