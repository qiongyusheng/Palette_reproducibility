library(tidyverse)
library(future)

FastSparseRowScale <- function(mat, scale = TRUE, center = TRUE, scale_max = 10, display_progress = TRUE) {
  .Call('_Seurat_FastSparseRowScale', PACKAGE = 'Seurat', mat, scale, center, scale_max, display_progress)
}

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

ScaleData1 <- function(
    object,
    features = NULL,
    do.scale = TRUE,
    do.center = TRUE,
    scale.max = 10,
    block.size = 1000,
    min.cells.to.block = 3000,
    verbose = TRUE,
    ...
) {
  features <- features %||% rownames(x = object)
  features <- as.vector(x = intersect(x = features, y = rownames(x = object)))
  object <- object[features, , drop = FALSE]
  object.names <- dimnames(x = object)
  min.cells.to.block <- min(min.cells.to.block, ncol(x = object))

  split.cells <- split(x = colnames(x = object), f = TRUE)
  gc(verbose = FALSE)
  
  if (verbose && (do.scale || do.center)) {
    msg <- paste(
      na.omit(object = c(
        ifelse(test = do.center, yes = 'centering', no = NA_character_),
        ifelse(test = do.scale, yes = 'scaling', no = NA_character_)
      )),
      collapse = ' and '
    )
    msg <- paste0(
      toupper(x = substr(x = msg, start = 1, stop = 1)),
      substr(x = msg, start = 2, stop = nchar(x = msg)),
      ' data matrix'
    )
    message(msg)
  }
  if (inherits(x = object, what = c('dgCMatrix', 'dgTMatrix'))) {
    scale.function <- FastSparseRowScale
  } else {
    object <- as.matrix(x = object)
    scale.function <- Seurat::FastRowScale
  }
  if (future::nbrOfWorkers() > 1) {
    blocks <- ChunkPoints1(dsize = length(x = features), csize = block.size)
    chunks <- expand.grid(
      names(x = split.cells),
      1:ncol(x = blocks),
      stringsAsFactors = FALSE
    )
    scaled.data <- future_lapply(
      X = 1:nrow(x = chunks),
      FUN = function(index) {
        row <- chunks[index, ]
        group <- row[[1]]
        block <- as.vector(x = blocks[, as.numeric(x = row[[2]])])
        arg.list <- list(
          mat = object[features[block[1]:block[2]], split.cells[[group]], drop = FALSE],
          scale = do.scale,
          center = do.center,
          scale_max = scale.max,
          display_progress = FALSE
        )
        arg.list <- arg.list[intersect(x = names(x = arg.list), y = names(x = formals(fun = scale.function)))]
        data.scale <- do.call(what = scale.function, args = arg.list)
        dimnames(x = data.scale) <- dimnames(x = object[features[block[1]:block[2]], split.cells[[group]]])
        suppressWarnings(expr = data.scale[is.na(x = data.scale)] <- 0)
        gc(verbose = FALSE)
        return(data.scale)
      }
    )
    if (length(x = split.cells) > 1) {
      merge.indices <- lapply(
        X = 1:length(x = split.cells),
        FUN = seq.int,
        to = length(x = scaled.data),
        by = length(x = split.cells)
      )
      scaled.data <- lapply(
        X = merge.indices,
        FUN = function(x) {
          return(suppressWarnings(expr = do.call(what = 'rbind', args = scaled.data[x])))
        }
      )
      scaled.data <- suppressWarnings(expr = do.call(what = 'cbind', args = scaled.data))
    } else {
      suppressWarnings(expr = scaled.data <- do.call(what = 'rbind', args = scaled.data))
    }
  } else {
    scaled.data <- matrix(
      data = NA_real_,
      nrow = nrow(x = object),
      ncol = ncol(x = object),
      dimnames = object.names
    )
    max.block <- ceiling(x = length(x = features) / block.size)
    for (x in names(x = split.cells)) {
      if (verbose) {
        if (length(x = split.cells) > 1 && (do.scale || do.center)) {
          message(gsub(pattern = 'matrix', replacement = 'from split ', x = msg), x)
        }
        pb <- txtProgressBar(min = 0, max = max.block, style = 3, file = stderr())
      }
      for (i in 1:max.block) {
        my.inds <- ((block.size * (i - 1)):(block.size * i - 1)) + 1
        my.inds <- my.inds[my.inds <= length(x = features)]
        arg.list <- list(
          mat = object[features[my.inds], split.cells[[x]], drop = FALSE],
          scale = do.scale,
          center = do.center,
          scale_max = scale.max,
          display_progress = FALSE
        )
        arg.list <- arg.list[intersect(x = names(x = arg.list), y = names(x = formals(fun = scale.function)))]
        data.scale <- do.call(what = scale.function, args = arg.list)
        dimnames(x = data.scale) <- dimnames(x = object[features[my.inds], split.cells[[x]]])
        scaled.data[features[my.inds], split.cells[[x]]] <- data.scale
        rm(data.scale)
        gc(verbose = FALSE)
        if (verbose) {
          setTxtProgressBar(pb = pb, value = i)
        }
      }
      if (verbose) {
        close(con = pb)
      }
    }
  }
  dimnames(x = scaled.data) <- object.names
  scaled.data[is.na(x = scaled.data)] <- 0
  gc(verbose = FALSE)
  return(scaled.data)
}


















