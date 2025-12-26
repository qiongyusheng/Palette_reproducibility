

FindClusters1 <- function(
    object,
    modularity.fxn = 1,
    initial.membership = NULL,
    node.sizes = NULL,
    resolution = 0.8,
    method = "matrix",
    algorithm = 1,
    n.start = 10,
    n.iter = 10,
    random.seed = 0,
    group.singletons = TRUE,
    temp.file.location = NULL,
    edge.file.name = NULL,
    verbose = TRUE,
    ...
) {
  #CheckDots(...)
  if (is.null(x = object)) {
    stop("Please provide an SNN graph")
  }
  if (tolower(x = algorithm) == "louvain") {
    algorithm <- 1
  }
  if (tolower(x = algorithm) == "leiden") {
    algorithm <- 4
  }
  if (nbrOfWorkers() > 1) {
    clustering.results <- future_lapply(
      X = resolution,
      FUN = function(r) {
        if (algorithm %in% c(1:3)) {
          ids <- RunModularityClustering1(
            SNN = object,
            modularity = modularity.fxn,
            resolution = r,
            algorithm = algorithm,
            n.start = n.start,
            n.iter = n.iter,
            random.seed = random.seed,
            print.output = verbose,
            temp.file.location = temp.file.location,
            edge.file.name = edge.file.name
          )
        } else if (algorithm == 4) {
          ids <- RunLeiden1(
            object = object,
            method = method,
            partition.type = "RBConfigurationVertexPartition",
            initial.membership = initial.membership,
            node.sizes = node.sizes,
            resolution.parameter = r,
            random.seed = random.seed,
            n.iter = n.iter
          )
        } else {
          stop("algorithm not recognised, please specify as an integer or string")
        }
        names(x = ids) <- colnames(x = object)
        ids <- GroupSingletons1(ids = ids, SNN = object, verbose = verbose)
        results <- list(factor(x = ids))
        names(x = results) <- paste0('res.', r)
        return(results)
      }
    )
    clustering.results <- as.data.frame(x = clustering.results)
  } else {
    clustering.results <- data.frame(row.names = colnames(x = object))
    for (r in resolution) {
      if (algorithm %in% c(1:3)) {
        ids <- RunModularityClustering1(
          SNN = object,
          modularity = modularity.fxn,
          resolution = r,
          algorithm = algorithm,
          n.start = n.start,
          n.iter = n.iter,
          random.seed = random.seed,
          print.output = verbose,
          temp.file.location = temp.file.location,
          edge.file.name = edge.file.name)
      } else if (algorithm == 4) {
        ids <- RunLeiden1(
          object = object,
          method = method,
          partition.type = "RBConfigurationVertexPartition",
          initial.membership = initial.membership,
          node.sizes = node.sizes,
          resolution.parameter = r,
          random.seed = random.seed,
          n.iter = n.iter
        )
      } else {
        stop("algorithm not recognised, please specify as an integer or string")
      }
      names(x = ids) <- colnames(x = object)
      ids <- GroupSingletons1(ids = ids, SNN = object, group.singletons = group.singletons, verbose = verbose)
      clustering.results[, paste0("res.", r)] <- factor(x = ids)
    }
  }
  return(clustering.results)
}



RunModularityClusteringCpp <- function(SNN, modularityFunction, resolution, algorithm, nRandomStarts, nIterations, randomSeed, printOutput, edgefilename) {
  .Call('_Seurat_RunModularityClusteringCpp', PACKAGE = 'Seurat', SNN, modularityFunction, resolution, algorithm, nRandomStarts, nIterations, randomSeed, printOutput, edgefilename)
}

RunModularityClustering1 <- function(
    SNN = matrix(),
    modularity = 1,
    resolution = 0.8,
    algorithm = 1,
    n.start = 10,
    n.iter = 10,
    random.seed = 0,
    print.output = TRUE,
    temp.file.location = NULL,
    edge.file.name = NULL
) {
  edge_file <- edge.file.name %||% ''
  clusters <- RunModularityClusteringCpp(
    SNN,
    modularity,
    resolution,
    algorithm,
    n.start,
    n.iter,
    random.seed,
    print.output,
    edge_file
  )
  return(clusters)
}


library(reticulate)
library(leiden)
RunLeiden1 <- function(
    object,
    method = c("matrix", "igraph"),
    partition.type = c(
      'RBConfigurationVertexPartition',
      'ModularityVertexPartition',
      'RBERVertexPartition',
      'CPMVertexPartition',
      'MutableVertexPartition',
      'SignificanceVertexPartition',
      'SurpriseVertexPartition'
    ),
    initial.membership = NULL,
    node.sizes = NULL,
    resolution.parameter = 1,
    random.seed = 0,
    n.iter = 10
) {
  if (!py_module_available(module = 'leidenalg')) {
    stop(
      "Cannot find Leiden algorithm, please install through pip (e.g. pip install leidenalg).",
      call. = FALSE
    )
  }
  switch(
    EXPR = method,
    "matrix" = {
      input <- as(object = object, Class = "matrix")
    },
    "igraph" = {
      input <- if (inherits(x = object, what = 'list')) {
        graph_from_adj_list(adjlist = object)
      } else if (inherits(x = object, what = c('dgCMatrix', 'matrix', 'Matrix'))) {
        if (inherits(x = object, what = 'Graph')) {
          object <- as.sparse(x = object)
        }
        graph_from_adjacency_matrix(adjmatrix = object, weighted = TRUE)
      } else if (inherits(x = object, what = 'igraph')) {
        object
      } else {
        stop(
          "Method for Leiden not found for class", class(x = object),
          call. = FALSE
        )
      }
    },
    stop("Method for Leiden must be either 'matrix' or igraph'")
  )
  #run leiden from CRAN package (calls python with reticulate)
  partition <- leiden(
    object = input,
    partition_type = partition.type,
    initial_membership = initial.membership,
    weights = NULL,
    node_sizes = node.sizes,
    resolution_parameter = resolution.parameter,
    seed = random.seed,
    n_iterations = n.iter
  )
  return(partition)
}




GroupSingletons1 <- function(ids, SNN, group.singletons = TRUE, verbose = TRUE) {
  # identify singletons
  singletons <- c()
  singletons <- names(x = which(x = table(ids) == 1))
  singletons <- intersect(x = unique(x = ids), singletons)
  if (!group.singletons) {
    ids[which(ids %in% singletons)] <- "singleton"
    return(ids)
  }
  # calculate connectivity of singletons to other clusters, add singleton
  # to cluster it is most connected to
  cluster_names <- as.character(x = unique(x = ids))
  cluster_names <- setdiff(x = cluster_names, y = singletons)
  connectivity <- vector(mode = "numeric", length = length(x = cluster_names))
  names(x = connectivity) <- cluster_names
  new.ids <- ids
  for (i in singletons) {
    i.cells <- names(which(ids == i))
    for (j in cluster_names) {
      j.cells <- names(which(ids == j))
      subSNN <- SNN[i.cells, j.cells]
      set.seed(1) # to match previous behavior, random seed being set in WhichCells
      if (is.object(x = subSNN)) {
        connectivity[j] <- sum(subSNN) / (nrow(x = subSNN) * ncol(x = subSNN))
      } else {
        connectivity[j] <- mean(x = subSNN)
      }
    }
    m <- max(connectivity, na.rm = T)
    mi <- which(x = connectivity == m, arr.ind = TRUE)
    closest_cluster <- sample(x = names(x = connectivity[mi]), 1)
    ids[i.cells] <- closest_cluster
  }
  if (length(x = singletons) > 0 && verbose) {
    message(paste(
      length(x = singletons),
      "singletons identified.",
      length(x = unique(x = ids)),
      "final clusters."
    ))
  }
  return(ids)
}
