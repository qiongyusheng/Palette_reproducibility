ComputeSNN <- function(nn_ranked, prune) {
  .Call('_Seurat_ComputeSNN', PACKAGE = 'Seurat', nn_ranked, prune)
}

FindNeighbors1 <- function(
    object,
    query = NULL,
    dist.matrix = FALSE, #是否返回距离矩阵
    k.param = 20,
    compute.SNN = TRUE,
    prune.SNN = 1/15,
    nn.method = "annoy",
    n.trees = 50,
    annoy.metric = "euclidean",
    nn.eps = 0,
    verbose = TRUE,
    force.recalc = FALSE,
    l2.norm = FALSE,
    cache.index = FALSE,
    index = NULL,
    ...
) {
  # CheckDots(...)
  if (is.null(x = dim(x = object))) {
    warning(
      "Object should have two dimensions, attempting to coerce to matrix",
      call. = FALSE
    )
    object <- as.matrix(x = object)
  }
  if (is.null(rownames(x = object))) {
    stop("Please provide rownames (cell names) with the input object")
  }
  n.cells <- nrow(x = object)
  if (n.cells < k.param) {
    warning(
      "k.param set larger than number of cells. Setting k.param to number of cells - 1.",
      call. = FALSE
    )
    k.param <- n.cells - 1
  }
  if (l2.norm) {
    object <- L2Norm1(mat = object)
    query <- query %iff% L2Norm1(mat = query)
  }
  query <- query %||% object
  # find the k-nearest neighbors for each single cell
  
  if (verbose) {
    message("Computing nearest neighbor graph")
  }
  nn.ranked <- NNHelper1(
    data = object,
    query = query,
    k = k.param,
    method = nn.method,
    n.trees = n.trees,
    searchtype = "standard",
    eps = nn.eps,
    metric = annoy.metric,
    cache.index = cache.index,
    index = index
  )
  nn.idx <- slot(object = nn.ranked, name = "nn.idx")
  rownames(x = nn.idx) <- slot(object = nn.ranked, name = "cell.names")
  
  # convert nn.ranked into a Graph
  j <- as.numeric(x = t(x = nn.idx))
  i <- ((1:length(x = j)) - 1) %/% k.param + 1
  nn.matrix <- as(sparseMatrix(i = i, j = j, 
                               x = 1, 
                               dims = c(nrow(x = query), 
                                        nrow(x = object))), Class = "Graph")
  rownames(x = nn.matrix) <- rownames(x = query)
  colnames(x = nn.matrix) <- rownames(x = object)
  neighbor.graphs <- list(nn = nn.matrix)
  
  if(dist.matrix){
    nn.dist <- slot(object = nn.ranked, name = "nn.dist")
    rownames(x = nn.dist) <- slot(object = nn.ranked, name = "cell.names")
    nnd.matrix <- as(sparseMatrix(i = i, j = j, 
                                  x = as.numeric(x = t(x = nn.dist)), 
                                  dims = c(nrow(x = query), 
                                           nrow(x = object))), Class = "Graph")
    # input <- as(object = object, Class = "matrix") Graph转换回matrix的方法
    rownames(x = nnd.matrix) <- rownames(x = query)
    colnames(x = nnd.matrix) <- rownames(x = object)
    neighbor.graphs[["nn.dist"]] <- nnd.matrix
  }
  
  if (compute.SNN) {
    if (verbose) {
      message("Computing SNN")
    }
    snn.matrix <- ComputeSNN(
      nn_ranked = nn.idx,
      prune = prune.SNN
    )
    rownames(x = snn.matrix) <- rownames(x = object)
    colnames(x = snn.matrix) <- rownames(x = object)
    snn.matrix <- SeuratObject::as.Graph(x = snn.matrix)
    neighbor.graphs[["snn"]] <- snn.matrix
  }
  return(neighbor.graphs)
}


NNHelper1 <- function(data, query = data, k, method, cache.index = FALSE, ...) {
  args <- as.list(x = sys.frame(which = sys.nframe()))
  args <- c(args, list(...))
  results <- (
    switch(
      EXPR = method,
      "rann" = {
        args <- args[intersect(x = names(x = args), y = names(x = formals(fun = nn2)))]
        do.call(what = 'nn2', args = args)
      },
      "annoy" = {
        args <- args[intersect(x = names(x = args), y = names(x = formals(fun = AnnoyNN1)))]
        do.call(what = 'AnnoyNN1', args = args)
      },
      stop("Invalid method. Please choose one of 'rann', 'annoy'")
    )
  )
  
  # 创建了一个定义的对象
  #############################################################
  n.ob <- new(
    Class = 'Neighbor',
    nn.idx = results$nn.idx,
    nn.dist = results$nn.dists,
    alg.info = results$alg.info %||% list(),
    cell.names = rownames(x = query)
  )
  #############################################################
  
  if (isTRUE(x = cache.index) && !is.null(x = results$idx)) {
    slot(object = n.ob, name = "alg.idx") <- results$idx
  }
  return(n.ob)
}


AnnoyNN1 <- function(data,
                    query = data,
                    metric = "euclidean",
                    n.trees = 50,
                    k,
                    search.k = -1,
                    include.distance = TRUE,
                    index = NULL
) {
  idx <- index %||% AnnoyBuildIndex1(
    data = data,
    metric = metric,
    n.trees = n.trees)
  nn <- AnnoySearch1(
    index = idx,
    query = query,
    k = k,
    search.k = search.k,
    include.distance = include.distance)
  nn$idx <- idx
  nn$alg.info <- list(metric = metric, ndim = ncol(x = data))
  return(nn)
}

library(RcppAnnoy)
AnnoyBuildIndex1 <- function(data, metric = "euclidean", n.trees = 50) {
  f <- ncol(x = data)
  a <- switch(
    EXPR = metric,
    "euclidean" =  new(Class = RcppAnnoy::AnnoyEuclidean, f),
    "cosine" = new(Class = RcppAnnoy::AnnoyAngular, f),
    "manhattan" = new(Class = RcppAnnoy::AnnoyManhattan, f),
    "hamming" = new(Class = RcppAnnoy::AnnoyHamming, f),
    stop ("Invalid metric")
  )
  for (ii in seq(nrow(x = data))) {
    a$addItem(ii - 1, data[ii, ])
  }
  a$build(n.trees)
  return(a)
}

library(future)
library(future.apply)
AnnoySearch1 <- function(index, query, k, search.k = -1, include.distance = TRUE) {
  n <- nrow(x = query)
  idx <- matrix(nrow = n,  ncol = k)
  dist <- matrix(nrow = n, ncol = k)
  convert <- methods::is(index, "Rcpp_AnnoyAngular")
  if (!inherits(x = plan(), what = "multicore")) {
    oplan <- plan(strategy = "sequential")
    on.exit(plan(oplan), add = TRUE)
  }
  res <- future_lapply(X = 1:n, FUN = function(x) {
    res <- index$getNNsByVectorList(query[x, ], k, search.k, include.distance)
    # Convert from Angular to Cosine distance
    if (convert) {
      res$dist <- 0.5 * (res$dist * res$dist)
    }
    list(res$item + 1, res$distance)
  })
  for (i in 1:n) {
    idx[i, ] <- res[[i]][[1]]
    if (include.distance) {
      dist[i, ] <- res[[i]][[2]]
    }
  }
  return(list(nn.idx = idx, nn.dists = dist))
}




L2Norm1 <- function(mat, MARGIN = 1){
  normalized <- Sweep1(
    x = mat,
    MARGIN = MARGIN,
    STATS = apply(
      X = mat,
      MARGIN = MARGIN,
      FUN = function(x){
        sqrt(x = sum(x ^ 2))
      }
    ),
    FUN = "/"
  )
  normalized[!is.finite(x = normalized)] <- 0
  return(normalized)
}
# Sweep out array summaries
#
# Reimplmentation of \code{\link[base]{sweep}} to maintain compatability with
# both R 3.X and 4.X
#
# @inheritParams base::sweep
# @param x an array.
#
# @seealso \code{\link[base]{sweep}}
#
Sweep1 <- function(x, MARGIN, STATS, FUN = '-', check.margin = TRUE, ...) {
  if (any(grepl(pattern = 'X', x = names(x = formals(fun = sweep))))) {
    return(sweep(
      X = x,
      MARGIN = MARGIN,
      STATS = STATS,
      FUN = FUN,
      check.margin = check.margin,
      ...
    ))
  } else {
    return(sweep(
      x = x,
      MARGIN = MARGIN,
      STATS = STATS,
      FUN = FUN,
      check.margin = check.margin,
      ...
    ))
  }
}