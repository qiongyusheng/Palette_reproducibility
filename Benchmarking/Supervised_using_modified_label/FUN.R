library(methods)
library(Matrix)
library(tidyverse)
setClassUnion(name = 'AnyMatrix', members = c("matrix", "dgCMatrix"))

Musa.Assay <- methods::setClass(
  Class = "Musa.Assay",
  slots = c(
    Modal = "character",
    Sample = "character",
    Raw = "AnyMatrix",
    data = "AnyMatrix",
    snn = "AnyMatrix",
    RepresentativeCell = "list",
    embedding = "AnyMatrix",
    var.features = "character",
    n.Count = "vector",
    n.Feature = "vector"
  )
)

Musa <- methods::setClass(
  Class = "Musa",
  slots = c(
    Assays = "list",
    Meta = "data.frame",
    Data.index = "list",
    Co.Features = "list",
    pre.clusters = "list",
    Subspaces = "list",
    emb.inf = "list",
    Int.result ="list" ,
    cluster = "character",
    Reduction = "list"
  )
)

Create.Musa.Object <- function(data.list, 
                               samples, 
                               modals,
                               filter = FALSE,
                               min.cells = 0, 
                               min.features = 0,
                               modal.order,
                               sparce = TRUE){
  if(missing(data.list) || missing(samples) || missing(modals)){
    stop("Must provide data.list, sample information and modal information.")
  }
  
  # 检查数据、批次、模态信息是否匹配
  if(length(data.list) != length(samples) || length(data.list) != length(modals) || length(samples) != length(modals)){
    stop("Wrong match information, data.list, samples and modals must same length")
  }
  
  #检查过滤信息
  if(filter){
    modal.order <- modal.order %||% unique(modals)
    a <- setdiff(modals,modal.order)
    if(length(a) != 0){
      stop("Incomplete modal information in modal.order")
    }
    b <- setdiff(modal.order,modals)
    if(length(b) != 0){
      stop("Redundant modal information in modal.order")
    }
    
    min.cells <- min.cells %||% as.numeric(rep(0,length(modal.order)))
    min.features <- min.features %||% as.numeric(rep(0,length(modal.order)))
    if(length(min.cells) != length(modal.order) || length(min.features) != length(modal.order)){
      stop("'min.cells' and 'min.features' must have same length as modal.order")
    }
  }
  
  # modals，sample和assay的名字都不能包含&和_
  if(Reduce("+",grepl('&', modals)) > 0){
    warning("character '&' in modal name will replaced by '-")
    modals <- gsub('&', '-', modals)
  }
  if(Reduce("+",grepl('_', modals)) > 0){
    warning("character '_' in modal name will replaced by '-")
    modals <- gsub('_', '-', modals)
  }
  if(Reduce("+",grepl('&', samples)) > 0){
    warning("character '&' in sample name will replaced by '-")
    samples <- gsub('&', '-', samples)
  }
  if(Reduce("+",grepl('_', samples)) > 0){
    warning("character '_' in sample name will replaced by '-")
    samples <- gsub('_', '-', samples)
  }
  
  # 检查数据有没有名字，没有的话用批次和模态信息命名
  if(is.null(names(data.list)) || any(names(x = data.list) == '')){
    names(data.list) <- paste0(modals,"_",samples)
  }else{
    if(Reduce("+",grepl('&', names(data.list))) > 0){
      warning("character '&' in data name will replaced by '-")
      names(data.list) <- gsub('&', '-', names(data.list))
    }
    if(Reduce("+",grepl('_', names(data.list))) > 0){
      warning("character '_' in data name will replaced by '-")
      names(data.list) <- gsub('_', '-', names(data.list))
    }
  }
  data.name <- names(data.list)
  
  # 检查每个数据矩阵
  # 每个数据矩阵是否存在，行名和列名是否存在，是否稀疏转换
  data.list <- lapply(data.list, check.data, sparce = sparce)
  
  # 模态内和批次内检查与过滤
  # 每个批次的细胞名称、数量、顺序是否对应
  # 每个模态的特征是否对应
  if(filter){
    min.cells <- min.cells[match(modals, modal.order)]
    min.features <- min.features[match(modals, modal.order)]
    
    filter.cell.list <- lapply(1:length(data.list), function(x){
      if(min.features[x] > 0){
        nfeatures <- Matrix::colSums(x = data.list[[x]] > 0)
        x <- colnames(data.list[[x]])[which(nfeatures >= min.features[x])]
      }else{
        x <- colnames(data.list[[x]])
      }
      x
    })
    # names(filter.cell.list) <- data.name
    filter.feature.list <- lapply(1:length(data.list), function(x){
      if(min.cells[x] >0){
        num.cells <- Matrix::rowSums(x = data.list[[x]] > 0)
        x <- rownames(data.list[[x]])[which(num.cells >= min.cells[x])]
      }else{
        x <- rownames(data.list[[x]])
      }
      x
    })
    # names(filter.feature.list) <- data.name
  }else{
    filter.cell.list <- lapply(data.list, colnames)
    filter.feature.list <- lapply(data.list, rownames)
  }
  
  # 在每个模态中保留共同特征
  for(j in unique(modals)){
    idx.modal <- which(modals == j)
    feature.list <- filter.feature.list[idx.modal]
    if(length(idx.modal) > 1){
      modal.feature <- Reduce(intersect, feature.list)
      data.list[idx.modal] <- lapply(1:length(idx.modal), function(x){
        x <- data.list[idx.modal][[x]][modal.feature, ]
        x
      })
    }
  }
  
  # 检查过滤之后细胞和特征是否存在没有了的情况
  dim.vector <- unlist(lapply(data.list,dim))
  if(is.numeric(Reduce("*",as.numeric(dim.vector))) && Reduce("*",as.numeric(dim.vector)) != 0){
    names(data.list) <- data.name
    assay <- lapply(1:length(data.list),function(x){
      x <- Create.Musa.Assay(count = data.list[[x]], 
                             modal = modals[x], 
                             sample = samples[x],
                             filter = FALSE,
                             sparce = sparce)
      x
    })
    names(assay) <- data.name
    data.idx <- list(assay.idx = data.frame(sample = samples, 
                                            modal = modals,class = rep("original",length(samples))))
    cell.name <- unlist(lapply(unique(samples), function(x){
      idx <- which(samples == x)[1]
      x <- colnames(data.list[[idx]])
      x
    }))
    Musa <- methods::new(
      Class = "Musa",
      Assays = assay,
      Meta = data.frame(cell.name),
      Data.index = data.idx
    )
  }else{
    stop("Check whether the filter conditions or the input information match")
  }
  return(Musa)
}

check.data <- function(counts, 
                       sparce = TRUE){
  
  if(missing(x = counts)){
    stop("Must provide 'counts'")
  }else if(!missing(x = counts)){
    if(anyDuplicated(x = rownames(x = counts))){
      warning(
        "Non-unique features (rownames) present in the input matrix, making unique",
        call. = FALSE,
        immediate. = TRUE
      )
      rownames(x = counts) <- make.unique(names = rownames(x = counts))
    }
    if(anyDuplicated(x = colnames(x = counts))){
      warning(
        "Non-unique cell names (colnames) present in the input matrix, making unique",
        call. = FALSE,
        immediate. = TRUE
      )
      colnames(x = counts) <- make.unique(names = colnames(x = counts))
    }
    if (is.null(x = colnames(x = counts))) {
      stop("No cell names (colnames) names present in the input matrix")
    }
    if (any(rownames(x = counts) == '')) {
      stop("Feature names of counts matrix cannot be empty", call. = FALSE)
    }
    if (nrow(x = counts) > 0 && is.null(x = rownames(x = counts))) {
      stop("No feature names (rownames) names present in the input matrix")
    }
    if(sparce){
      if (!inherits(x = counts, what = 'dgCMatrix')){
        counts <- as(as.matrix(counts), 'dgCMatrix')
      }
    }
  }
  return(counts)
}

Create.Musa.Assay <- function(count, 
                              filter = TRUE,
                              min.cells = 0, 
                              min.features = 0, 
                              modal = "RNA",
                              sample = "Batch_1",
                              slot = "Raw",
                              sparce = TRUE){
  count <- check.data(counts = count, sparce = sparce)
  if(filter){
    # filter genes on the number of cells expressing
    if (min.cells > 0) {
      num.cells <- Matrix::rowSums(x = count > 0)
      count <- count[which(x = num.cells >= min.cells), ]
    }
    # Filter based on min.features
    if (min.features > 0) {
      nfeatures <- Matrix::colSums(x = count > 0)
      count <- count[, which(x = nfeatures >= min.features)]
    }
  }
  # Ensure row- and column-names are vectors, not arrays
  if (!is.vector(x = rownames(x = count))) {
    rownames(x = count) <- as.vector(x = rownames(x = count))
  }
  if (!is.vector(x = colnames(x = count))) {
    colnames(x = count) <- as.vector(x = colnames(x = count))
  }
  n.Count <- Matrix::colSums(x = count)
  n.Feature <- Matrix::colSums(x = count > 0)
  if(slot == "Raw"){
    assay <- methods::new(
      Class = "Musa.Assay",
      Modal = modal,
      Sample = sample,
      Raw = count,
      n.Count = n.Count,
      n.Feature = n.Feature
    )
  }else if(slot == "data"){
    assay <- methods::new(
      Class = "Musa.Assay",
      Modal = modal,
      Sample = sample,
      Raw = count,
      data = count,
      n.Count = n.Count,
      n.Feature = n.Feature
    )
  }else{
    stop("Unknown slot.")
  }
  
  return(assay)
}

Add.Assay <- function(assay, musa.object, check = TRUE, class = "original",name = NULL){
  modal <- assay@Modal
  sample <- assay@Sample
  modals <- musa.object@Data.index[["assay.idx"]][,"modal"]
  samples <- musa.object@Data.index[["assay.idx"]][,"sample"]
  classes <- musa.object@Data.index[["assay.idx"]][,"class"]
  
  modal_logic <- modals == modal
  sample_logic <- samples == sample
  
  if(Reduce("+",modal_logic) == 0 && Reduce("+",sample_logic) == 0){
    warning("No modal or sample match between Musa object and new assay")
    new.name <- colnames(assay@Raw)
    if(dim(musa.object@Meta)[2] > 1){
      musa.object@Meta <- musa.object@Meta[,"cell.name"]
      warning("Length mismatched meta data removed")
    }
    cell.name <- c(musa.object@Meta[,"cell.name"],new.name)
    musa.object@Meta <- data.frame(cell.name)
  }else{
    if(check){
      if(Reduce("+",modal_logic) != 0){
        idx.m <- which(modal_logic == 1)
        feature.list <- list(old.feature = rownames(musa.object@Assays[[idx.m[1]]]@Raw),
                             new.feature = rownames(assay@Raw))
        co.feature <- Reduce(intersect, feature.list)
        for(i in idx.m){
          musa.object@Assays[[i]]@Raw <- musa.object@Assays[[i]]@Raw[co.feature,]
        }
        assay@Raw <- assay@Raw[co.feature,]
      }
      if(Reduce("+",sample_logic) != 0){
        idx.s <- which(sample_logic == 1)
        cell.list <- list(old.cell = colnames(musa.object@Assays[[idx.s[1]]]@Raw),
                          new.cell = colnames(assay@Raw))
        co.cell <- Reduce(intersect, cell.list)
        for(i in idx.s){
          musa.object@Assays[[i]]@Raw <- musa.object@Assays[[i]]@Raw[ , co.cell]
        }
        assay@Raw <- assay@Raw[,co.cell]
      }
    }
    # musa.object@Assays <- c(musa.object@Assays, assay)
  }
  name <- name %||%paste0(modal,"_",sample)
  musa.object@Assays[[name]] <- assay
  modals <- c(modals,modal)
  samples <- c(samples,sample)
  classes <- c(classes,class)
  data.idx <- data.frame(sample = samples, 
                         modal = modals,
                         class = classes)
  musa.object@Data.index[["assay.idx"]] <- data.idx
  return(musa.object)
}

############## Palette integration code ##############
library(dplyr)
library(patchwork)
library(tidyverse)
library(RSpectra)
library(doParallel)
library(foreach)
library(Seurat)
library(SeuratObject)
library(Rcpp)
library(uwot)
library(bigstatsr)


Normalize.Data <- function(object, 
                           modal = "RNA", 
                           normal.method = "LogNormalize", 
                           scale.factor = 10000, 
                           margin = 1){
  # x: an object
  # modal : which modality used for normalize
  # normal.method: RNA("LogNormalize"), ADT("CLR")
  # scale.factor: Sets the scale factor for cell-level normalization
  # margin: If performing CLR normalization, normalize across features (1) or cells (2)
  if(!(inherits(modal,what = "character"))){
    stop("Parameter modal should input as character vector.")
  }
  # !(inherits(x = dims, what = 'list'))
  if(!(inherits(normal.method,what = "character"))){
    stop("Parameter normal.method should input as character vector.")
  }
  
  if(!(inherits(margin,what = "numeric"))){
    stop("Parameter margin should input as numeric vector.")
  }
  
  #检查数据是否一一对应
  if(length(modal)!=length(normal.method) || length(normal.method) != length(margin)){
    stop("Unmatched information.")
  }
  if(!all(modal %in% object@Data.index[["assay.idx"]]$modal)){
    stop("Unknown modal information.")
  }
  if(!all(normal.method %in% c("LogNormalize","CLR"))){
    stop("Unknown normalize method.")
  }
  if(!all(margin %in% c(1,2))){
    stop("The margin parameter must be selected between 1 and 2.")
  }
  for(i in 1:length(modal)){
    modal_logic <- ifelse(object@Data.index[["assay.idx"]]$modal == modal[i], TRUE, FALSE)
    idx <- which(modal_logic)
    if(normal.method[i]=="CLR"){
      for(j in idx){
        object@Assays[[j]]@data <- as(NormalizeData1(object@Assays[[j]]@Raw, 
                                                     normalization.method = normal.method[i],
                                                     scale.factor = scale.factor,
                                                     margin = margin[i],
                                                     block.size = NULL,
                                                     verbose = FALSE),"dgCMatrix")
      }
    }else{
      for(j in idx){
        object@Assays[[j]]@data <- NormalizeData1(object@Assays[[j]]@Raw, 
                                                  normalization.method = normal.method[i],
                                                  scale.factor = scale.factor,
                                                  margin = margin[i],
                                                  block.size = NULL,
                                                  verbose = FALSE)
      }
    }
  }
  return(object)
}


Find.HVFs <- function(object, 
                      modal,
                      loess.span = 0.3,
                      clip.max = 'auto',
                      nfeatures = 2000,
                      quantile_norm = FALSE,
                      verbose = TRUE){
  
  idx <- which(object@Data.index[["assay.idx"]]$modal == modal)
  HVF.set <- c()
  for(i in idx){
    HVFs <- FindVariableFeatures1(object@Assays[[i]]@Raw,
                                  loess.span = loess.span,
                                  clip.max = clip.max,
                                  nfeatures = nfeatures,
                                  verbose = verbose)
    HVF.set <- c(HVF.set, HVFs)
    rm(HVFs)
  }
  co_feature <- unique(HVF.set)
  object@Co.Features[[modal]] <- co_feature
  
  if(quantile_norm == TRUE){
    name.list <- lapply(idx, function(x){
      x <- colnames(object@Assays[[x]]@Raw)
      x
    })
    data <- Reduce(cbind,lapply(idx,function(x){
      x <- object@Assays[[x]]@data[co_feature,]
      x
    }))
    q_data <- preprocessCore::normalize.quantiles(data)
    for(i in idx){
      object@Assays[[i]]@data <- q_data[,name.list[[i]]]
    }
  }
  return(object)
}

Get.Features <- function(object = NULL,
                         modal = NULL,
                         sample = NULL,
                         check = TRUE){
  if(is.null(object) || is.null(modal)){
    stop("Must provide Musa object and modal information")
  }else{
    if(is.null(sample)){
      modal_logic <- ifelse(object@Data.index[["assay.idx"]]$modal == modal, TRUE, FALSE)
      if(check){
        feature.set <- c()
        for(i in 1:length(modal_logic)){
          if(modal_logic[i]){
            feature <- rownames(object@Assays[[i]]@Raw)
            feature.set <- c(feature.set, feature)
            rm(feature)
          }
        }
        out <- unique(feature.set)
      }else{
        idx <- which(modal_logic == 1)[1]
        out <- rownames(object@Assays[[idx]]@Raw)
      }
    }else{
      modal_logic <- ifelse(object@Data.index[["assay.idx"]]$modal == modal, TRUE, FALSE)
      sample_logic <- ifelse(object@Data.index[["assay.idx"]]$sample == sample, TRUE, FALSE)
      idx <- which(modal_logic*sample_logic == 1)
      if(length(idx) != 1){
        stop("Please check the input modal and sample information")
      }else{
        out <- rownames(object@Assays[[idx]]@Raw)
      }
    }
  }
  return(out)
}

Add.HVFs <- function(object,
                     modal = NULL,
                     features = NULL,
                     quantile_norm = FALSE){
  if(missing(object) || is.null(modal)){
    stop("Must specify an object and provide modal information")
  }
  # modal_logic <- ifelse(object@Data.index[["assay.idx"]]$modal == modal, TRUE, FALSE)
  available.feautres <- Get.Features(object = object,
                                     modal = modal,
                                     check = TRUE)
  if(is.null(features)){
    message("no features provide, use all features by default")
    input.features <- available.feautres
  }else{
    input.features <- intersect(features, available.feautres)
    if(length(input.features)==0){
      stop("Provided features are not included in this modal")
    }
  }
  object@Co.Features[[modal]] <- input.features
  
  if(quantile_norm == TRUE){
    idx <- which(object@Data.index[["assay.idx"]]$modal == modal)
    name.list <- lapply(idx, function(x){
      x <- colnames(object@Assays[[x]]@Raw)
      x
    })
    data <- Reduce(cbind,lapply(idx,function(x){
      x <- object@Assays[[x]]@data[input.features,]
      x
    }))
    q_data <- preprocessCore::normalize.quantiles(data)
    for(i in idx){
      object@Assays[[i]]@data <- q_data[,name.list[[i]]]
    }
  }
  
  return(object)
}


Binary.Data <- function(object,
                        modal = NULL){
  
  if(missing(object) || missing(modal)){
    stop("Must provide a Musa object and modal information")
  }
  if(!(modal %in% object@Data.index[["assay.idx"]]$modal)){
    stop("Unknown modal information.")
  }
  for(i in 1:length(modal)){
    
    modali <- modal[i]
    modal_logic <- object@Data.index[["assay.idx"]]$modal %in% modali
    for(j in 1:length(modal_logic)){
      
      if(modal_logic[j]){
        data <- object@Assays[[j]]@Raw
        
        if (inherits(x = data, what = "CsparseMatrix")) {
          slot(object = data, name = "x") <- rep.int(
            x = 1,
            times = length(
              x = slot(object = data, name = "x")
            )
          )
          
        } else {
          data[data > 1] <- 1
        }
      }
    }
  }
  object@Assays[[j]]@data <- data
  return(object)
}

TFIDF.Data <- function(object,
                       modal,
                       sample = NULL,
                       scale.factor = 1e4,
                       idf = NULL,
                       verbose = TRUE){
  
  if(missing(object) || missing(modal) || length(modal) != 1){
    stop("Must provide a Musa object and a modal information")
  }
  if(!(modal %in% object@Data.index[["assay.idx"]]$modal)){
    stop("Unknown modal information.")
  }
  modal_logic <- object@Data.index[["assay.idx"]]$modal %in% modal
  
  if(!is.null(sample) && !(sample %in% object@Data.index[["assay.idx"]]$sample)){
    stop("Unknown samople information.")
  }
  sample <- sample %||% unique(object@Data.index[["assay.idx"]]$sample)
  sample_logic <- object@Data.index[["assay.idx"]]$sample %in% sample
  logic <- modal_logic * sample_logic
  
  if (verbose) {
    message("Performing TF-IDF normalization")
  }
  
  if(is.null(idf)){
    rsums <- list()
    col_sum <- 0
    for(i in 1:length(logic)){
      if(logic[i]){
        rsums <- c(rsums,list(rowSums(x = object@Assays[[i]]@Raw)))
        col_sum <- col_sum + ncol(object@Assays[[i]]@Raw)
      }
    }
    rsums.vec <- Reduce("+",rsums)
    idf <- col_sum / rsums.vec
  }else{
    idx <- whcih(logic)[1]
    if(length(idf) != nrow(x = object@Assays[[idx]]@Raw)){
      stop("Length of supplied IDF vector does not match",
           " number of rows in input matrix")
    }
  }
  
  for(j in 1:length(logic)){
    if(logic[j]){
      data <- object@Assays[[j]]@Raw
      if (inherits(x = data, what = "data.frame")) {
        data <- as.matrix(x = data)
      }
      if (!inherits(x = data, what = "CsparseMatrix")) {
        data <- as(object = data, Class = "CsparseMatrix")
      }
      npeaks <- colSums(x = data)
      if (any(npeaks == 0)) {
        warning("Some cells contain 0 total counts")
      }
      tf <- tcrossprod(x = data, y = Matrix::Diagonal(x = 1 / npeaks))
      norm.data <- Matrix::Diagonal(n = length(x = idf), x = idf)  %*% tf# ！！！！
      slot(object = norm.data, name = "x") <- log1p(
        x = slot(object = norm.data, name = "x") * scale.factor)
      # norm.data <- t(L2Norm1(t(norm.data)))
      colnames(x = norm.data) <- colnames(x = data)
      rownames(x = norm.data) <- rownames(x = data)
      # set NA values to 0
      vals <- slot(object = norm.data, name = "x")
      vals[is.na(x = vals)] <- 0
      slot(object = norm.data, name = "x") <- vals
      object@Assays[[j]]@data <- norm.data
    }
  }
  return(object)
}

IDFLog_try <- function(object,
                        modal,
                        sample = NULL,
                        scale.factor = 1e4,
                        idf = NULL,
                        verbose = TRUE){
  
  if(missing(object) || missing(modal) || length(modal) != 1){
    stop("Must provide a Musa object and a modal information")
  }
  if(!(modal %in% object@Data.index[["assay.idx"]]$modal)){
    stop("Unknown modal information.")
  }
  modal_logic <- object@Data.index[["assay.idx"]]$modal %in% modal
  
  if(!is.null(sample) && !(sample %in% object@Data.index[["assay.idx"]]$sample)){
    stop("Unknown samople information.")
  }
  sample <- sample %||% unique(object@Data.index[["assay.idx"]]$sample)
  sample_logic <- object@Data.index[["assay.idx"]]$sample %in% sample
  logic <- modal_logic * sample_logic
  
  if (verbose) {
    message("Performing IDF and log normalization")
  }
  
  if(is.null(idf)){
    rsums <- list()
    col_sum <- 0
    for(i in 1:length(logic)){
      if(logic[i]){
        rsums <- c(rsums,list(rowSums(x = object@Assays[[i]]@Raw)))
        col_sum <- col_sum + ncol(object@Assays[[i]]@Raw)
      }
    }
    rsums.vec <- Reduce("+",rsums)
    idf <- col_sum / rsums.vec
  }else{
    idx <- whcih(logic)[1]
    if(length(idf) != nrow(x = object@Assays[[idx]]@Raw)){
      stop("Length of supplied IDF vector does not match",
           " number of rows in input matrix")
    }
  }
  
  for(j in 1:length(logic)){
    if(logic[j]){
      data <- object@Assays[[j]]@Raw
      if (inherits(x = data, what = "data.frame")) {
        data <- as.matrix(x = data)
      }
      if (!inherits(x = data, what = "CsparseMatrix")) {
        data <- as(object = data, Class = "CsparseMatrix")
      }
      norm.data <- Matrix::Diagonal(n = length(x = idf), x = idf)%*%data
      # slot(object = norm.data, name = "x") <- log1p(
      #   x = slot(object = norm.data, name = "x") * scale.factor)
      norm.data <- NormalizeData1(norm.data,
                                  normalization.method = "LogNormalize",
                                  scale.factor = scale.factor,
                                  margin = 1,
                                  block.size = NULL,
                                  verbose = FALSE)
      colnames(x = norm.data) <- colnames(x = data)
      rownames(x = norm.data) <- rownames(x = data)
      # set NA values to 0
      vals <- slot(object = norm.data, name = "x")
      vals[is.na(x = vals)] <- 0
      slot(object = norm.data, name = "x") <- vals
      object@Assays[[j]]@data <- norm.data
    }
  }
  return(object)
}

LogIDFL2_try <- function(object,
                          modal,
                          sample = NULL,
                          scale.factor = 1e4,
                          idf = NULL,
                          verbose = TRUE){
  
  if(missing(object) || missing(modal) || length(modal) != 1){
    stop("Must provide a Musa object and a modal information")
  }
  if(!(modal %in% object@Data.index[["assay.idx"]]$modal)){
    stop("Unknown modal information.")
  }
  modal_logic <- object@Data.index[["assay.idx"]]$modal %in% modal
  
  if(!is.null(sample) && !(sample %in% object@Data.index[["assay.idx"]]$sample)){
    stop("Unknown samople information.")
  }
  sample <- sample %||% unique(object@Data.index[["assay.idx"]]$sample)
  sample_logic <- object@Data.index[["assay.idx"]]$sample %in% sample
  logic <- modal_logic * sample_logic
  
  if (verbose) {
    message("Performing LogIDF and L2 normalization")
  }
  sample <- sample %||% unique(object@Data.index[["assay.idx"]]$sample)
  sample_logic <- object@Data.index[["assay.idx"]]$sample %in% sample
  logic <- modal_logic * sample_logic
  
  if (verbose) {
    message("Performing IDF and log normalization")
  }
  
  if(is.null(idf)){
    rsums <- list()
    col_sum <- 0
    for(i in 1:1:length(logic)){
      if(logic[i]){
        rsums <- c(rsums,list(rowSums(x = object@Assays[[i]]@Raw)))
        col_sum <- col_sum + ncol(object@Assays[[i]]@Raw)
      }
    }
    rsums.vec <- Reduce("+",rsums)
    idf <- log1p(col_sum / rsums.vec)
  }else{
    idx <- whcih(logic)[1]
    if(length(idf) != nrow(x = object@Assays[[idx]]@Raw)){
      stop("Length of supplied IDF vector does not match",
           " number of rows in input matrix")
    }
  }
  
  for(j in 1:length(logic)){
    if(logic[j]){
      data <- object@Assays[[j]]@Raw
      if (inherits(x = data, what = "data.frame")) {
        data <- as.matrix(x = data)
      }
      if (!inherits(x = data, what = "CsparseMatrix")) {
        data <- as(object = data, Class = "CsparseMatrix")
      }
      norm.data <- Matrix::Diagonal(n = length(x = idf), x = idf) %*% data
      L2.vec <- sparse_L2(norm.data@p,norm.data@x)
      norm.data@x <- as.vector(L2.vec)
      colnames(x = norm.data) <- colnames(x = data)
      rownames(x = norm.data) <- rownames(x = data)
      # set NA values to 0
      vals <- slot(object = norm.data, name = "x")
      vals[is.na(x = vals)] <- 0
      slot(object = norm.data, name = "x") <- vals
      object@Assays[[j]]@data <- norm.data
    }
  }
  return(object)
  
}


Run_Multi_SVD <- function(x, n.dims = 50, thresh = 10000L){
  
  if(nrow(x[[1]]) < thresh){
    
    d <- lapply(x, function(a){
      a <- tcrossprod(a)/(ncol(a)-1)
      a
    })
    
    d <- Reduce("+", d)
    V <- irlba::irlba(d, nv = n.dims, work = min(1000,n.dims*3))[["u"]] %>% as.matrix()
    r <- lapply(x, function(b){
      crossprod(V, b)
    })
    return(r)
    
  }else{
    
    if(Reduce("+",lapply(x,ncol)) > nrow(x[[1]])){
      d <- bigstatsr::FBM(nrow(x[[1]]),nrow(x[[1]]), init = 0)
      for(i in 1:length(x)){
        tmpx <- x[[i]]
        tmp <- bigstatsr::FBM(nrow(tmpx),ncol(tmpx),init = 0)
        ind_nozero <- which(tmpx != 0, arr.ind = TRUE)
        vec <- tmpx@x/sqrt(ncol(tmpx)-1)
        tmpx@x <- as.vector(vec)
        tmp[ind_nozero] <- tmpx[ind_nozero]
        d <- bigstatsr::big_increment(d,bigstatsr::big_tcrossprodSelf(tmp)[])
      }
      
      V <- bigstatsr::big_SVD(d,k = n.dims)[["u"]] %>% as.matrix()
      r <- lapply(x, function(b){
        crossprod(V, b)
      })
      return(r)
        
      }else{
      d <- bigstatsr::FBM(nrow(x[[1]]),Reduce("+",lapply(x,ncol)), init = 0)
      cells <- c(0,cumsum(lapply(x,ncol)))
      for(i in 1:length(x)){
        tmpx <- x[[i]]
        vec <- tmpx@x/sqrt(ncol(tmpx)-1)
        tmpx@x <- as.vector(vec)
        ind_nozero <- which(tmpx != 0, arr.ind = TRUE)
        ind2 <- ind_nozero
        ind2[,2] <- ind2[,2]+cells[i]
        d[ind2] <- tmpx[ind_nozero]
      }
    
      V <- bigstatsr::big_SVD(d,k = n.dims)[["u"]] %>% as.matrix()
      r <- lapply(x, function(b){
        crossprod(V, b)
      })
      return(r)
    
      }
    }
}


Run_Multi_SVD_s <- function(x, n.dims = 50, thresh = 10000L){
  
  if(nrow(x[[1]]) < thresh){
    
    d <- lapply(x, function(a){
      a <- tcrossprod(a)/(ncol(a)-1)
      a
    })
    
    d <- Reduce("+", d)
    V <- RSpectra::svds(d, k = n.dims)[["u"]] %>% as.matrix()
    r <- lapply(x, function(b){
      crossprod(V, b)
    })
    return(list(r,V))
    
  }else{
    d <- lapply(x, function(a){
      a <- a/sqrt(ncol(a)-1)
      a
    })
    d <- Reduce(cbind, d)
    V <- RSpectra::svds(d, k = n.dims)[["u"]] %>% as.matrix()
    r <- lapply(x, function(b){
      crossprod(V, b)
    })
    return(list(r,V))
  }
}

FastSparseRowScaleWithKnownStats <- function(mat, mu, sigma, scale = TRUE, center = TRUE, scale_max = 10, display_progress = TRUE) {
  .Call('_Seurat_FastSparseRowScaleWithKnownStats', PACKAGE = 'Seurat', mat, mu, sigma, scale, center, scale_max, display_progress)
}

Find.Cluster <- function(object,
                         modal = NULL,
                         sample = NULL,
                         dim.reduce = TRUE,
                         dims = NULL,
                         method = c("PCA","LSI"), 
                         joint = FALSE,
                         resolution = 1,
                         nn.k = 20,
                         nn.method = "annoy",
                         annoy.metric = "euclidean",
                         prune.SNN = 1/15,
                         downsample = TRUE,
                         quantile.cutoff = 0.75,
                         min.cell = 20,
                         verbose = TRUE){
  if(missing(object)){
    stop("Must provide an object information")
  }
  
  modal <- modal %||% unique(object@Data.index[["assay.idx"]]$modal)
  sample <- sample %||% unique(object@Data.index[["assay.idx"]]$sample)
  
  if(length(resolution) != length(modal)){
    resolution <- rep(resolution[1],length(modal))
  }
  
  if(dim.reduce){
    
    if(length(setdiff(method,c("PCA","LSI")))){
      stop("Unknown dimensionality reduction method.")
    }
    
    if(!(inherits(x = dims, what = 'list'))){
      dims <- list(dims)
    }
    
    if(length(dims) != length(modal)){
      dims <- lapply(1:length(modal),function(x){
        x <- dims[[1]]
      })
      # names(dims) <- modal
    }
    names(dims) <- modal
    
    if(length(method) == 1){
      method <- rep(method,length(modal))
      # names(method) <- modal
    }else if(length(method) != length(modal)){
      stop("Unmatched length between parameters method and modal.")
    }
    names(method) <- modal
  }
  
  modal_logic <- object@Data.index[["assay.idx"]]$modal %in% modal
  sample_logic <- object@Data.index[["assay.idx"]]$sample %in% sample
  logic <- modal_logic * sample_logic
  
  if(Reduce("+",logic) == 0){
    stop("Musa Assay does not exist")
  }
  
  cluster.name <- names(object@Assays)
  for(i in 1:length(modal)){
    
    idx_i <- which(object@Data.index[["assay.idx"]]$modal == modal[i])
    cluster.name_i <- cluster.name[idx_i]
    features <- object@Co.Features[[modal[i]]]
    if(is.null(features)){
      stop("请先鉴定高可变异特征！！")
    }
    data.list_i <- lapply(idx_i, function(x) object@Assays[[x]]@data[features,])
    names(data.list_i) <- cluster.name_i
    if(dim.reduce){
      
      method_i <- method[i]
      dim_i <- dims[[i]]
      
      if(method_i == "PCA"){
        scale.list <- lapply(data.list_i,function(x){
          x <- ScaleData1(object = x,
                          features = features,
                          do.scale = TRUE,
                          do.center = TRUE,
                          scale.max = 10,
                          block.size = 1000,
                          min.cells.to.block = 3000,
                          verbose = FALSE)
          x
        })
        names(scale.list) <- cluster.name_i
        if(joint){
          input.data <- Run_Multi_SVD(scale.list,
                                      max(dim_i))
        }else{
          input.data <- lapply(1:length(scale.list),function(x){
            V <- irlba::irlba(scale.list[[x]], nv = max(dim_i), work = max(dim_i)*3)[["u"]] %>% as.matrix()
            x <- crossprod(V, scale.list[[x]])
            x
          })
        }
        input.data <- lapply(input.data, function(x) x[dim_i,])
        
        names(input.data) <- cluster.name_i
        rm(data.list_i,scale.list)
        gc(verbose = FALSE)
        
      }else if(method_i == "LSI"){
        
        if(joint){
          input.data <- Run_Multi_SVD(data.list_i,
                                      max(dim_i))
        }else{
          
          if(nrow(data.list_i[[1]]) > 10000L){
            input.data <- lapply(1:length(data.list_i),function(x){
              d <- bigstatsr::FBM(nrow(data.list_i[[x]]),
                                  ncol(data.list_i[[x]]),
                                  init = 0)
              ind_nozero <- which(data.list_i[[x]] != 0, arr.ind = TRUE)
              d[ind_nozero] <- data.list_i[[x]][ind_nozero]
              V <- bigstatsr::big_SVD(d,k = max(dim_i))[["u"]] %>% as.matrix()
              x <- crossprod(V, data.list_i[[x]])
              x
            })
            
          }else{
            input.data <- lapply(1:length(data.list_i),function(x){
              V <- irlba::irlba(data.list_i[[x]], nv = max(dim_i), work = max(dim_i)*3)[["u"]] %>% as.matrix()
              x <- crossprod(V, data.list_i[[x]])
              x
            })
          }
          
        }
        input.data <- lapply(input.data, function(x) x[dim_i,])
        names(input.data) <- cluster.name_i
        rm(data.list_i)
        gc(verbose = FALSE)
      }
      
    }else{
      input.data <- data.list_i
      rm(data.list_i)
      gc(verbose = FALSE)
    }
    
    out_i <- list()
    for(j in 1:length(input.data)){
      snn <- FindNeighbors1(t(input.data[[j]]),
                            k.param = nn.k,
                            nn.method = nn.method,
                            prune.SNN = prune.SNN,
                            annoy.metric = annoy.metric,
                            verbose = verbose)$snn
      
      out_i[[j]] <- FindClusters1(snn,
                                  resolution = resolution[i],
                                  verbose = verbose)
      
      Cluster <- as.character(out_i[[j]][,1])
      
      if(downsample == TRUE){
        ASW_score <- silhouette_cpp(as.numeric(out_i[[j]][,1]),as.matrix(input.data[[j]]))
        Cluster.uni <- as.character(unique(Cluster))
        flag <- c()
        if(quantile.cutoff > 0){
          for(jj in 1:length(Cluster.uni)){
            Index <- which(Cluster == Cluster.uni[jj])
            asw_jj <- ASW_score[Index]
            if(length(asw_jj) < min.cell){
             flag <- c(flag,Index) 
            }else{
              thresh <- stats::quantile(asw_jj,quantile.cutoff)
              index <- which(asw_jj > thresh)
              if (length(index) > min.cell) {
                flag <- c(flag, Index[index])
              }else{
                index <- rank(asw_jj,ties.method = "max")[seq(min.cell)]
                flag <- c(flag, Index[index])
              }
            }
          }
        }else{
          flag <- seq(ncol(input.data[[j]]))
        }
        
        object@Assays[[cluster.name_i[j]]]@RepresentativeCell[["CellName"]] <- colnames(input.data[[j]])[flag]
        object@Assays[[cluster.name_i[j]]]@RepresentativeCell[["Cluster"]] <- Cluster[flag]
        
      }else{
        flag <- seq(ncol(input.data[[j]]))
        
        object@Assays[[cluster.name_i[j]]]@RepresentativeCell[["CellName"]] <- colnames(input.data[[j]])[flag]
        object@Assays[[cluster.name_i[j]]]@RepresentativeCell[["Cluster"]] <- Cluster[flag]
      }
    }
    names(out_i) <- cluster.name_i
    
    current.cluster <- names(object@pre.clusters)
    if(!is.null(current.cluster) && length(current.cluster) > 0){
      replace.idx <- which(current.cluster %in% cluster.name_i)
      if(length(replace.idx) > 0){
        object@pre.clusters <- object@pre.clusters[-replace.idx]
      }
    }
    object@pre.clusters <- c(object@pre.clusters,out_i)
  }
  
  return(object)
}

clust_labels_to_int = function(clust_labels) {
  uniq_labels = unique(clust_labels)
  clust_labels_int = rep(NA, length(clust_labels))
  
  for (ii in 1:length(clust_labels)) {
    clust_labels_int[ii] = match(clust_labels[ii], uniq_labels)
  }
  
  return(clust_labels_int)
}

DownSample <- function(object,
                       modal = NULL,
                       sample = NULL,
                       supervised = TRUE,
                       group = NULL,
                       dim.reduce = TRUE,
                       dims = NULL,
                       method = c("PCA","LSI"), # adt和rna选PCA，atac选LSI
                       joint = FALSE,
                       quantile.cutoff = 0.75,
                       min.cell = 20){
  
  # supervised: 是否是有监督的，如果是有监督的，则需要提供‘group’，其是位于object@Meta中的列名；
  #                             如果不是有监督的，则对应数据的object@pre.clusters应该已经包含了聚类信息
  
  
  # 使用Meta信息进行采样，适用于有标签情景
  if(missing(object)){
    stop("Must provide an object information")
  }
  
  if(supervised && is.null(group)){
    stop('Must selcet meta data under supervised case.')
  }
  if(supervised){
    if(length(intersect(group,colnames(object@Meta))) == 0){
      stop('Group information does not exist.')
    }
  }
  
  modal <- modal %||% unique(object@Data.index[["assay.idx"]]$modal)
  sample <- sample %||% unique(object@Data.index[["assay.idx"]]$sample)
  
  if(dim.reduce){
    
    if(length(setdiff(method,c("PCA","LSI")))){
      stop("Unknown dimensionality reduction method.")
    }
    
    if(!(inherits(x = dims, what = 'list'))){
      dims <- list(dims)
    }
    
    if(length(dims) != length(modal)){
      dims <- lapply(1:length(modal),function(x){
        x <- dims[[1]]
      })
      # names(dims) <- modal
    }
    names(dims) <- modal
    
    if(length(method) == 1){
      method <- rep(method,length(modal))
      # names(method) <- modal
    }else if(length(method) != length(modal)){
      stop("Unmatched length between parameters method and modal.")
    }
    names(method) <- modal
  }
  
  modal_logic <- object@Data.index[["assay.idx"]]$modal %in% modal
  sample_logic <- object@Data.index[["assay.idx"]]$sample %in% sample
  logic <- modal_logic * sample_logic
  
  if(Reduce("+",logic) == 0){
    stop("Musa Assay does not exist")
  }
  
  
  if(supervised){
    meta_cluster <- as.character(clust_labels_to_int(as.character(object@Meta[,group])))
    
    ####### 缺失细胞类型代码
    idx <- which(!is.na(meta_cluster))
    meta_cluster <- meta_cluster[idx]
    meta_cluster <- factor(meta_cluster,
                           levels = as.character(1:length(unique(meta_cluster))))
    names(meta_cluster) <- object@Meta[,'cell.name'][idx]
    #######
    
    # 以前的代码
    # meta_cluster <- factor(meta_cluster,
    #                        levels = as.character(1:length(unique(meta_cluster))))
    # names(meta_cluster) <- object@Meta[,'cell.name']
    
    for(i in 1:length(modal)){
      idx_i <- which(object@Data.index[["assay.idx"]]$modal == modal[i])
      for(j in idx_i){
        clust_name <- names(object@Assays)[j]
        tmp_name <- colnames(object@Assays[[j]]@Raw)
        
        ####### 缺失细胞类型代码
        tmp_name <- intersect(tmp_name,names(meta_cluster))
        #######
        
        tmp_clust <- as.data.frame(meta_cluster[tmp_name])
        rownames(tmp_clust) <- tmp_name
        object@pre.clusters[[clust_name]] <- tmp_clust
      }
    }
  }else{
    for(i in 1:length(modal)){
      idx_i <- which(object@Data.index[["assay.idx"]]$modal == modal[i])
      for(j in idx_i){
        clust_name <- paste0(modal[i],'_',object@Data.index[["assay.idx"]]$sample[j])
        if(is.null(object@pre.clusters[[clust_name]])){
          stop('Clustering information must be provided')
        }
      }
    }
  }
  
  cluster.name <- names(object@Assays)
  for(i in 1:length(modal)){
    
    idx_i <- which(object@Data.index[["assay.idx"]]$modal == modal[i])
    cluster.name_i <- cluster.name[idx_i]
    features <- object@Co.Features[[modal[i]]]
    if(is.null(features)){
      stop("请先鉴定高可变异特征！！")
    }
    data.list_i <- lapply(idx_i, function(x) object@Assays[[x]]@data[features,])
    
    idx <- which(names(object@pre.clusters) %in% cluster.name_i)
    idx <- which(names(object@pre.clusters) %in% cluster.name_i)
    Cluster_list <- lapply(idx,function(x){
      x <- object@pre.clusters[[x]]
      x
    })
    names(Cluster_list) <- names(object@pre.clusters)[idx]
    Cluster_list <- Cluster_list[cluster.name_i]
    
    names(data.list_i) <- cluster.name_i
    if(dim.reduce){
      
      method_i <- method[i]
      dim_i <- dims[[i]]
      
      if(method_i == "PCA"){
        scale.list <- lapply(data.list_i,function(x){
          x <- ScaleData1(object = x,
                          features = features,
                          do.scale = TRUE,
                          do.center = TRUE,
                          scale.max = 10,
                          block.size = 1000,
                          min.cells.to.block = 3000,
                          verbose = FALSE)
          x
        })
        names(scale.list) <- cluster.name_i
        if(joint){
          input.data <- Run_Multi_SVD(scale.list,
                                      max(dim_i))
        }else{
          input.data <- lapply(1:length(scale.list),function(x){
            V <- irlba::irlba(scale.list[[x]], nv = max(dim_i), work = max(dim_i)*3)[["u"]] %>% as.matrix()
            x <- crossprod(V, scale.list[[x]])
            x
          })
        }
        input.data <- lapply(input.data, function(x) x[dim_i,])
        
        names(input.data) <- cluster.name_i
        rm(data.list_i,scale.list)
        gc(verbose = FALSE)
        
      }else if(method_i == "LSI"){
        
        if(joint){
          input.data <- Run_Multi_SVD(data.list_i,
                                      max(dim_i))
        }else{
          
          input.data <- lapply(1:length(data.list_i),function(x){
            V <- RSpectra::svds(data.list_i[[x]], k = max(dim_i))[["u"]] %>% as.matrix()
            # irlba::irlba(data.list_i[[x]], nv = max(dim_i), work = max(dim_i)*3)[["u"]] %>% as.matrix()
            x <- crossprod(V, data.list_i[[x]])
            x
          })
          
        }
        input.data <- lapply(input.data, function(x) x[dim_i,])
        names(input.data) <- cluster.name_i
        rm(data.list_i)
        gc(verbose = FALSE)
      }
      
    }else{
      input.data <- data.list_i
      rm(data.list_i)
      gc(verbose = FALSE)
    }
    
    for(j in 1:length(input.data)){
      
      Cluster <- as.character(Cluster_list[[j]][,1])
      ASW_score <- silhouette_cpp(as.numeric(Cluster_list[[j]][,1]),as.matrix(input.data[[j]]))
      Cluster.uni <- as.character(unique(Cluster))
      flag <- c()
      if(quantile.cutoff > 0){
        for(jj in 1:length(Cluster.uni)){
          Index <- which(Cluster == Cluster.uni[jj])
          asw_jj <- ASW_score[Index]
          if(length(asw_jj) < min.cell){
            flag <- c(flag,Index)
          }else{
            thresh <- stats::quantile(asw_jj,quantile.cutoff)
            index <- which(asw_jj > thresh)
            if (length(index) > min.cell) {
              flag <- c(flag, Index[index])
            }else{
              index <- rank(asw_jj,ties.method = "max")[seq(min.cell)]
              flag <- c(flag, Index[index])
            }
          }
        }
      }else{
        flag <- seq(ncol(input.data[[j]]))
      }
      
      object@Assays[[cluster.name_i[j]]]@RepresentativeCell[["CellName"]] <- colnames(input.data[[j]])[flag]
      object@Assays[[cluster.name_i[j]]]@RepresentativeCell[["Cluster"]] <- Cluster[flag]
    }
    current.cluster <- names(object@pre.clusters)
    if(!is.null(current.cluster) && length(current.cluster) > 0){
      replace.idx <- which(current.cluster %in% cluster.name_i)
      if(length(replace.idx) > 0){
        object@pre.clusters <- object@pre.clusters[-replace.idx]
      }
    }
    object@pre.clusters <- c(object@pre.clusters,Cluster_list)
  }
  
  return(object)
  
}


scale_data <- function(data.x, 
                       do.center = T, 
                       do.scale = T, 
                       row.means = NULL, 
                       row.sds = NULL) {
  if (do.center) {
    if (is.null(row.means)) {
      data_mean <- Matrix::rowMeans(data.x)
    } else {
      data_mean <- row.means
    }
    data.x <- data.x - sapply(1:ncol(data.x), function(i) data_mean)
  }
  if (do.scale) {
    if (is.null(row.sds)) {
      data_stddev <- matrixStats::rowSds(as.matrix(data.x))
    } else {
      data_stddev <- row.sds
    }
    index <- which(data_stddev > 0)
    data.x[index, ] <- data.x[index, ] / sapply(1:ncol(data.x), function(i) data_stddev[index])
  }
  data.x
}

Build.Kernel.cluster3 <- function(data.list,
                                  cluster,
                                  # snn,
                                  dim.reduce = TRUE,
                                  cos.dims = NULL,
                                  Angle.var = 10,
                                  max.Angle = 35,
                                  joint = TRUE,
                                  seed = 123L){
  
  set.seed(seed)
  
  cluster.uni <- list()
  idx <- list()
  for(i in 1:length(cluster)){
    cluster_i <- cluster[[i]][,1]
    clusteri.uni <- as.character(unique(cluster_i))
    ni <- length(clusteri.uni)
    idxi <- lapply(1:ni,function(x) which(cluster_i == clusteri.uni[x]))
    cluster.uni[[i]] <- clusteri.uni
    idx[[i]] <- idxi
  }
  
  Angle.var <- Angle.var*pi/180
  max.Angle <- max.Angle*pi/180
  
  
  message("Performing batch specific cluster median cosine similarity. \n")
  gc(verbose = F)
  
  if(joint && dim.reduce){
    
    subKernel <- list()
    scaler <- c()
    for(i in 1:(length(cluster)-1)){
      for(j in (i+1):length(cluster)){
        tmp.list <- data.list[c(i,j)]
        svd.list <- Run_Multi_SVD(tmp.list, max(cos.dims))
        svd.list <- lapply(svd.list, function(x) as.matrix(x[cos.dims,]))
        L2.list <- lapply(svd.list, function(x) L2Norm1(t(x)))
        rm(svd.list)
        
        tmpi <- FindMedian_self(tcrossprod(L2.list[[1]]),idx[[i]])
        rownames(tmpi) <- colnames(tmpi) <- cluster.uni[[i]]
        
        tmpj <- FindMedian_self(tcrossprod(L2.list[[2]]),idx[[j]])
        rownames(tmpj) <- colnames(tmpj) <- cluster.uni[[j]]
        gc(verbose = F)
        
        tmp <- FindMedian(L2.list[[1]]%*%t(L2.list[[2]]),idx[[i]],idx[[j]])
        filterLink <- FindGroupLink(tmp,
                                    tmpi,
                                    tmpj,
                                    Angle.var,
                                    max.Angle)
        
        link <- PairwiseKernel_norm(L2.list[[1]]%*%t(L2.list[[2]]),
                                    filterLink,
                                    idx[[i]],
                                    idx[[j]],
                                    seed = seed)
        
        subKernel <- c(subKernel,link[['subK']])
        scaler <- c(scaler,link[["scaler"]])
        
      }
    }
    
  }else{
    if(dim.reduce){
      message("Performing SVD and L2 norm. \n")
      svd.list <- Run_Multi_SVD(data.list, max(cos.dims))
      svd.list <- lapply(svd.list, function(x) as.matrix(x[cos.dims,]))
      L2.list <- lapply(svd.list, function(x) L2Norm1(t(x)))
      message("SVD and L2 norm done! \n")
      rm(svd.list)  
    }else{
      message("Performing L2 norm. \n")
      L2.list <- lapply(data.list, function(x) L2Norm1(t(x)))
      message("L2 norm done! \n")
    }
    gc(verbose = F)
    
    self.sim <- list()
    
    message("Performing batch specific cluster median cosine similarity. \n")
    for(i in 1:length(cluster)){
      tmp <- FindMedian_self(tcrossprod(L2.list[[i]]),idx[[i]])
      rownames(tmp) <- colnames(tmp) <- cluster.uni[[i]]
      self.sim[[names(data.list)[i]]] <- tmp
    }
    gc(verbose = F)
    
    subKernel <- list()
    scaler <- c()
    for(i in 1:(length(cluster)-1)){
      for(j in (i+1):length(cluster)){
        tmp <- FindMedian(L2.list[[i]]%*%t(L2.list[[j]]),idx[[i]],idx[[j]])
        filterLink <- FindGroupLink(tmp,
                                    self.sim[[i]],
                                    self.sim[[j]],
                                    Angle.var,
                                    max.Angle)
        link <- PairwiseKernel_norm(L2.list[[i]]%*%t(L2.list[[j]]),
                                    filterLink,
                                    idx[[i]],
                                    idx[[j]],
                                    seed = seed)
        
        subKernel <- c(subKernel,link[['subK']])
        scaler <- c(scaler,link[["scaler"]])
      }
    }
    
  }
  
  gc(verbose = F)
  max_s <- max(scaler)
  
  N <- Reduce("+",lapply(data.list,ncol))
  ZK <- Matrix::Matrix(0, nrow = N, ncol = N)

  K <- Insert_submat(ZK,subKernel,sapply(data.list,ncol),scale = T,scaler = max_s)
  K <- Matrix::drop0(K)
  
  Name.list <- lapply(data.list,colnames)
  TName <- Reduce("c",Name.list)
  colnames(K) <- rownames(K) <- TName

  gc(verbose = F)
  
  return(K)
}

Build.Kernel.cluster3_S <- function(data.list,
                                    cluster,
                                    seed = 123L){
  
  set.seed(seed)
  
  cluster.uni <- unique(unlist(cluster))
  nn <- length(cluster.uni)
  idx <- list()
  for(i in 1:length(cluster)){
    cluster_i <- cluster[[i]][,1]
    idxi <- lapply(1:nn,function(x){
      tmp_idx <- which(cluster_i == cluster.uni[x])
      if(length(tmp_idx) == 0){
        tmp_idx <- -1
      }
      x <- tmp_idx
      x
    })
    idx[[i]] <- idxi
  }

  N_list <- sapply(data.list,ncol)
  N <- Reduce("+",N_list)
  K_list  <- sKernel_norm(N_list,
                     N,nn,
                     idx)
  K <- K_list[['K']]
  B <- K_list[['B']]
  rm(K_list)
  K <- Matrix::drop0(K)
  B <- Matrix::drop0(B)
  
  Name.list <- lapply(data.list,colnames)
  TName <- Reduce("c",Name.list)
  colnames(K) <- rownames(K) <- TName
  colnames(B) <- rownames(B) <- TName
  
  return(list(K,B))
  
}


Find_SubSpace3 <- function(object,
                          modal = NULL,
                          sample = NULL,
                          lambda_i = NULL,
                          dim.reduce = TRUE,
                          sub_dims = NULL,
                          cos_dims = NULL,
                          Angle_var = 10,
                          max_Angle = 35,
                          joint = TRUE,
                          proj_scale = T,
                          proj_center = T,
                          L2.norm = TRUE,
                          seed = 123L,
                          pre.reduce = TRUE,
                          pre_dims = list(1:300),
                          RP_thresh = 10000,
                          ...){
  
  set.seed(seed)
  
  idx <- which(object@Data.index[["assay.idx"]]$modal == modal)
  features <- object@Co.Features[[modal]]
  data.input <- lapply(idx,function(x) object@Assays[[x]]@data[features,object@Assays[[x]]@RepresentativeCell[[1]]])
  cluster.input <- lapply(idx,function(x) as.data.frame(object@Assays[[x]]@RepresentativeCell[[2]]))
  names(data.input) <- names(object@Assays)[idx]
  
  K <- Build.Kernel.cluster3(data.input,
                                  cluster.input,
                                  joint = joint,
                                  dim.reduce = TRUE,
                                  cos.dims = cos_dims,
                                  Angle.var = Angle_var,
                                  max.Angle = max_Angle,
                                  seed = seed)
  
  data.mat <- Reduce(cbind,lapply(idx,function(x) object@Assays[[x]]@data[features,]))
  if(nrow(data.mat) > RP_thresh){
    if(pre.reduce){
      message(paste0("Performing pre-dimensionality reduction on modal ",modal,". \n"))
      ncells <- sum(unlist(lapply(idx,function(x) ncol(object@Assays[[x]]@data))))
      pre_dims[[1]] <- as.integer(intersect(1:ncells,pre_dims[[1]]))
      data_res <- Run_Multi_SVD(lapply(idx,function(x) object@Assays[[x]]@data[features,]),
                                n.dims = max(pre_dims[[1]]))
      
      data_res <- as.matrix(Reduce(cbind,data_res))[as.integer(pre_dims[[1]]),]
      colnames(data_res) <- colnames(data.mat)
      rownames(data_res) <- paste0("dim_",seq(nrow(data_res)))
      rm(data.mat)
      gc(verbose = F)
    }else{
      message(paste0("Performing random projection on modal ",modal,". \n"))
      RP_mat <- JL_Trans(nrow(data.mat),seed = seed)
      data_FBM <- bigstatsr::FBM(nrow(data.mat), ncol(data.mat), init = 0)
      ind_nozero <- which(data.mat != 0, arr.ind = TRUE)
      data_FBM[ind_nozero] <- data.mat[ind_nozero]
      # data_res <- bigstatsr::big_cprodMat(RP_mat,data_FBM)[]
      data_res <- bigstatsr::big_apply(data_FBM,function(X,ind){
        big_cprodMat(RP_mat, X[, ind])
      },a.combine = "cbind",block.size = 1000)[]
      
      colnames(data_res) <- colnames(data.mat)
      rownames(data_res) <- paste0("dim_",seq(nrow(data_res)))
      rm(data.mat)
      gc(verbose = F)
    }
  }else{
    data_res <- data.mat
    rm(data.mat)
    gc(verbose = F)
  }
  scaled.data <- ScaleData1(object = data_res,
                            features = rownames(data_res),
                            do.scale = proj_scale,
                            do.center = proj_center,
                            scale.max = 10,
                            block.size = 1000,
                            min.cells.to.block = 3000,
                            verbose = FALSE)
  
  data.use <- scaled.data[,colnames(K)]
  DVec <- rowSums(K)
  fill <- mean(DVec[which(DVec != 0)])
  Zidx <- which(diag(K) == 0)
  diag(K)[Zidx] <- fill*(1-lambda_i)
  
  r <- GetS(data.use, K, DVec, lambda_i)
  Subspace <- suppressWarnings(RSpectra::eigs_sym(calAv_cpp,k=max(sub_dims),which = "LA",n = nrow(r),args = r))[["vectors"]] %>% as.matrix() # 使用函数，保证数值稳定
  L_dim <- as.integer(intersect(seq(ncol(Subspace)),sub_dims))
  Subspace <- Subspace[,L_dim]
  embedding <- crossprod(as.matrix(scaled.data), Subspace)
  if(exists('RP_mat')){
    Subspace <- RP_mat%*%Subspace
    rm(RP_mat)
  }

  if(L2.norm){
    embedding <- L2Norm1(embedding)
  }
  rownames(embedding) <- colnames(scaled.data)
  colnames(embedding) <- colnames(Subspace) <- paste0("dim_",seq(ncol(embedding)))
  
  return(list(emb = embedding, Index = idx, Subspace = Subspace))
}

Find_SubSpace3_s <- function(object,
                             modal = NULL,
                             sample = NULL,
                             lambda_i = NULL,
                             dim.reduce = TRUE,
                             sub_dims = NULL,
                             joint = TRUE,
                             proj_scale = T,
                             proj_center = T,
                             L2.norm = TRUE,
                             seed = 123L,
                             pre.reduce = TRUE,
                             pre_dims = list(1:300),
                             RP_thresh = 10000,
                             ...){
  
  set.seed(seed)
  
  idx <- which(object@Data.index[["assay.idx"]]$modal == modal)
  features <- object@Co.Features[[modal]]
  # data.input <- lapply(idx,function(x) object@Assays[[x]]@data[features,])
  data.input <- lapply(idx,function(x) object@Assays[[x]]@data[features,object@Assays[[x]]@RepresentativeCell[[1]]])
  cluster.input <- lapply(idx,function(x) as.data.frame(object@Assays[[x]]@RepresentativeCell[[2]]))
  names(data.input) <- names(object@Assays)[idx]
  
  
  K_list <- Build.Kernel.cluster3_S(data.input,
                             cluster.input,
                             seed = seed)
  
  data.mat <- Reduce(cbind,lapply(idx,function(x) object@Assays[[x]]@data[features,]))
  if(nrow(data.mat) > RP_thresh){
    if(pre.reduce){
      message(paste0("Performing pre-dimensionality reduction on modal ",modal,". \n"))
      ncells <- sum(unlist(lapply(idx,function(x) ncol(object@Assays[[x]]@data))))
      pre_dims[[1]] <- as.integer(intersect(1:ncells,pre_dims[[1]]))
      res <- Run_Multi_SVD_s(lapply(idx,function(x) object@Assays[[x]]@data[features,]),
                                n.dims = max(pre_dims[[1]]))
      data_res <- res[[1]]
      V_pre <- res[[2]]
      rm(res)
      rownames(V_pre) <- features
      data_res <- as.matrix(Reduce(cbind,data_res))[as.integer(pre_dims[[1]]),]
      V_pre <- V_pre[,as.integer(pre_dims[[1]])]
      colnames(data_res) <- colnames(data.mat)
      rownames(data_res) <- paste0("dim_",seq(nrow(data_res)))
      colnames(V_pre) <- paste0("dim_",seq(nrow(data_res)))
      rm(data.mat)
      gc(verbose = F)
    }else{
      # 执行随机投影
      message(paste0("Performing random projection on modal ",modal,". \n"))
      RP_mat <- JL_Trans(nrow(data.mat),seed = seed)
      data_FBM <- bigstatsr::FBM(nrow(data.mat), ncol(data.mat), init = 0)
      ind_nozero <- which(data.mat != 0, arr.ind = TRUE)
      data_FBM[ind_nozero] <- data.mat[ind_nozero]
      # data_res <- bigstatsr::big_cprodMat(RP_mat,data_FBM)[]
      data_res <- bigstatsr::big_apply(data_FBM,function(X,ind){
        big_cprodMat(RP_mat, X[, ind])
      },a.combine = "cbind",block.size = 1000)[]
      
      colnames(data_res) <- colnames(data.mat)
      rownames(data_res) <- paste0("dim_",seq(nrow(data_res)))
      rm(data.mat)
      gc(verbose = F)
    }
  }else{
    data_res <- data.mat
    rm(data.mat)
    gc(verbose = F)
  }
  
  scaled.data <- ScaleData1(object = data_res,
                            features = rownames(data_res),
                            do.scale = proj_scale,
                            do.center = proj_center,
                            scale.max = 10,
                            block.size = 1000,
                            min.cells.to.block = 3000,
                            verbose = FALSE)
  
  K <- K_list[[1]]
  B <- K_list[[2]]
  rm(K_list)
  data.use <- scaled.data[,colnames(K)]
  
  r <- GetS_new(data.use, K, B, lambda_i)
  Subspace <- suppressWarnings(RSpectra::eigs_sym(calAv_cpp,k=max(sub_dims),which = "LA",n = nrow(r),args = r))[["vectors"]] %>% as.matrix() # 使用函数，保证数值稳定
  L_dim <- as.integer(intersect(seq(ncol(Subspace)),sub_dims)) # 如果数值不稳定取不到规定维数的子空间，则用得到的子空间
  Subspace <- Subspace[,L_dim]
  embedding <- crossprod(as.matrix(scaled.data), Subspace)
  # if(exists('RP_mat')){
  #   Subspace <- RP_mat%*%Subspace
  #   rm(RP_mat)
  # }
  # if(exists('V')){
  #   Subspace <- V%*%Subspace
  #   rm(RP_mat)
  # }
  
  if(L2.norm){
    embedding <- L2Norm1(embedding)
  }
  rownames(embedding) <- colnames(scaled.data)
  colnames(embedding) <- colnames(Subspace) <- paste0("dim_",seq(ncol(embedding)))
  # rownames(Subspace) <- features
  
  return(list(emb = embedding, Index = idx, Subspace = Subspace))
}

Find.Subspace3 <- function(object,
                          modal = NULL,
                          sample = NULL,
                          lambda = c(0.9),#调节批次内结构保留程度的参数，取值在0到1之间
                          # alpha = c(0.2),# 调节细胞间的距离，使原始数据上没有SNN连接的细胞保持较远距离
                          dim.reduce = TRUE,
                          sub.dims = NULL,
                          cos.dims = NULL,
                          Angle.var = 10,
                          max.Angle = 35,
                          joint = TRUE,
                          supervised = FALSE,
                          proj.scale = T,
                          proj.center = T,
                          L2.norm = TRUE,
                          seed = 123L,
                          pre.reduce = FALSE,
                          pre_dims = list(1:3000),
                          RP_thresh = 10000,
                          ...){
  set.seed(seed)
  
  if(missing(object)){
    stop("Must provide object")
  }
  sample <- sample %||% unique(object@Data.index[["assay.idx"]]$sample)
  modal <- modal %||% unique(object@Data.index[["assay.idx"]]$modal)
  
  # 去掉只包含单一批次的数据
  drop.idx <- c()
  for(i in 1:length(modal)){
    tmp_idx <- which(object@Data.index[["assay.idx"]]$modal == modal[i])
    if(length(tmp_idx) == 1){
      drop.idx <- c(drop.idx,i)
    }
  }
  
  if(length(drop.idx) > 0){
    modal <- modal[-drop.idx]
  }
  
  if(!(inherits(sub.dims,what = "list"))){
    sub.dims <- list(sub.dims)
  }
  if(!(inherits(cos.dims,what = "list"))){
    cos.dims <- list(cos.dims)
  }
  
  if(length(dim.reduce) != length(modal)){
    dim.reduce <- rep(dim.reduce[1],length(modal))
  }
  if(length(sub.dims) != length(modal)){
    sub.dims <- replicate(length(modal),sub.dims[[1]], simplify = FALSE)
  }
  if(length(cos.dims) != length(modal)){
    cos.dims <- replicate(length(modal),cos.dims[[1]], simplify = FALSE)
  }
  if(!(inherits(lambda,what = "numeric"))){
    lambda <- as.numeric(lambda)
  }
  if(length(lambda) != length(modal)){
    lambda <- rep(lambda[1],length(modal))
  }
  
  if(supervised){
    
    for(i in 1:length(modal)){
      
      out <- Find_SubSpace3_s(object = object,
                            modal = modal[i],
                            sample = sample,
                            lambda_i = lambda[i],
                            dim.reduce = dim.reduce[i],
                            sub_dims = sub.dims[[i]],
                            joint = joint,
                            proj_scale = proj.scale,
                            proj_center = proj.center,
                            L2.norm = L2.norm,
                            seed = seed,
                            pre.reduce = pre.reduce,
                            pre_dims = pre_dims,
                            RP_thresh = RP_thresh,
                            ...)
      
      embedding <- out[["emb"]]
      idx <- out[["Index"]]
      Subspace <- out[["Subspace"]]
      
      for(j in 1:length(idx)){
        object@Assays[[idx[j]]]@embedding <- t(embedding[colnames(object@Assays[[idx[j]]]@data),])
      }
      object@Subspaces[[modal[i]]] <- Subspace
    }
    
  }else{
    # if(!(inherits(alpha,what = "numeric"))){
    #    alpha <- as.numeric(alpha)
    # }
    if(!(inherits(Angle.var,what = "numeric"))){
      Angle.var <- as.numeric(Angle.var)
    }
    if(!(inherits(max.Angle,what = "numeric"))){
      max.Angle <- as.numeric(max.Angle)
    }
    
    # if(!(inherits(euc.dims,what = "list"))){
    #   euc.dims <- list(euc.dims)
    # }
    # if(length(alpha) != length(modal)){
    #   alpha <- rep(alpha[1],length(modal))
    # }
    if(length(Angle.var) != length(modal)){
      Angle.var <- rep(Angle.var[1],length(modal))
    }
    if(length(max.Angle) != length(modal)){
      max.Angle <- rep(max.Angle[1],length(modal))
    }
    
    for(i in 1:length(modal)){
      
      out <- Find_SubSpace3(object = object,
                            modal = modal[i],
                            sample = sample,
                            lambda_i = lambda[i],
                            # alpha_i = alpha[i],
                            dim.reduce = dim.reduce[i],
                            sub_dims = sub.dims[[i]],
                            cos_dims = cos.dims[[i]],
                            Angle_var = Angle.var[i],
                            max_Angle = max.Angle[i],
                            joint = joint,
                            proj_scale = proj.scale,
                            proj_center = proj.center,
                            L2.norm = L2.norm,
                            seed = seed,
                            pre.reduce = pre.reduce,
                            pre_dims = pre_dims,
                            RP_thresh = RP_thresh,
                            ...)
      
      embedding <- out[["emb"]]
      idx <- out[["Index"]]
      Subspace <- out[["Subspace"]]
      
      for(j in 1:length(idx)){
        object@Assays[[idx[j]]]@embedding <- t(embedding[colnames(object@Assays[[idx[j]]]@data),])
      }
      object@Subspaces[[modal[i]]] <- Subspace
    }
  }
  
  return(object)
}

JL_Trans <- function(Dim,
                     eps = 0.1,
                     max_dim = 15000L,
                     seed = NA_integer_, 
                     method = "gaussian"){
  
  # method should be one of 'li' or 'gaussian'
  
  if (!is.na(x = seed)) {
    set.seed(seed = seed)
  }
  method <- method[1L]
  
  if (!is.null(x = eps)) {
    if (eps > 1 || eps <= 0) {
      stop("'eps' must be 0 < eps <= 1")
    }
    ncol <- floor(x = 4 * log(x = Dim) / ((eps ^ 2) / 2 - (eps ^ 3 / 3)))
  }
  
  ncols <- min(ncol, max_dim)
  
  if(ncols > Dim){
    message("Using raw dimensions.")
    random_matrix <- bigstatsr::FBM(Dim,Dim,init = 0)
    diag(random_matrix) <- 1
    return(random_matrix)
  }else{
    
    if(!(method %in% c("li","gaussian"))){
      stop("Only methods 'li' and 'gaussian' are supported.")
    }
    
    if(method == "li"){
      s <- ceiling(sqrt(ncols))
      pro <- c((1/(2*s)),(1-(1/s)),(1/(2*s)))
      # random_matrix <- matrix(sample(seq.int(-1L, 1L),size = Dim * ncols,replace = TRUE,prob = prob),nrow = nrow)
      random_matrix <- bigstatsr::FBM(Dim,ncols,init = sample(seq.int(-1L, 1L),size = Dim * ncols,replace = TRUE,prob = prob))
      return(random_matrix)
    }
    
    if(method == "gaussian"){
      #crandom_matrix <- matrix(rnorm(Dim*ncols,mean=0,sd=(1/sqrt(ncols))),Dim,ncols)
      random_matrix <- bigstatsr::FBM(Dim,ncols,init = rnorm(Dim*ncols,mean=0,sd=(1/sqrt(ncols))))
      return(random_matrix)
    }
  }
  
}

dataLink <- function(object){
  
  modal <- object@Data.index[["assay.idx"]]$modal
  sample <- object@Data.index[["assay.idx"]]$sample
  modal.uni <- unique(modal)
  sample.uni <- unique(sample)
  adj_mat <- as(matrix(0,length(modal.uni)+length(sample.uni),length(modal.uni)+length(sample.uni)),"dgCMatrix")
  rownames(adj_mat) <- colnames(adj_mat) <- c(modal.uni,sample.uni)
  
  for(i in 1:length(modal)){
    adj_mat[modal[i],sample[i]] <- adj_mat[sample[i],modal[i]] <- 1
  }
  
  g <- graph.adjacency(adj_mat, mode="undirected")
  if (components(g)$no != 1) {
    stop("modal network is not connected")
  }
  return(g)
}

IsolatedModalEmbedding <- function(object,
                                   modal = NULL,
                                   method = "PCA",
                                   dims = NULL,
                                   do.scale = TRUE){
  
 modal <- modal %||% check_one(object@Data.index[["assay.idx"]]$modal)
 modal <- check_aga(modal,object@Data.index[["assay.idx"]]$modal)
 if(length(modal) == 0){
   message("There are no isolated modal data.")
   return(object)
 }else{
   if(!(inherits(method,what = "character"))){
     method <- as.character(method)
   }
   if(!(inherits(dims,what = "list"))){
     dims <- list(dims)
   }
   if(length(method) != length(modal)){
     method <- rep(method[1],length(modal))
   }
   if(length(dims) != length(modal)){
     dims <- replicate(length(modal),dims[[1]], simplify = FALSE)
   }
   if(length(do.scale) != length(modal)){
     do.scale <- rep(do.scale[1],length(modal))
   }
   
   for(i in 1:length(modal)){
     idx <- which(object@Data.index[["assay.idx"]]$modal == modal[i])
     features <- object@Co.Features[[modal[i]]]
     data <- object@Assays[[idx]]@data[features,]
     if(method[i] == "PCA"){
       
       input <- ScaleData1(object = data,
                           features = features,
                           do.scale = do.scale[i],
                           do.center = TRUE,
                           scale.max = 10,
                           block.size = 1000,
                           min.cells.to.block = 3000,
                           verbose = FALSE)
       
     }else if(method[i] == "SVD"){
       
       if(do.scale[i]){
         input <- ScaleData1(object = data,
                             features = features,
                             do.scale = do.scale[i],
                             do.center = FALSE,
                             scale.max = 10,
                             block.size = 1000,
                             min.cells.to.block = 3000,
                             verbose = FALSE)
       }else{
         input <- data
       }
       
     }else{
       stop("Unknown dim reduce method.")
     }
     rm(data)
     V <- irlba::irlba(input, nv = max(dims[[i]]), work = max(dims[[i]])*3)[["u"]] %>% as.matrix()
     V <- V[,dims[[i]]]
     emb <- t(V) %*% as.matrix(input)
     colnames(emb) <- colnames(input)
     rownames(emb) <- colnames(V) <- paste0("dim_",seq(nrow(emb)))
     rownames(V) <- rownames(input)
     
     object@Assays[[idx]]@embedding <- emb
     object@Subspaces[[modal[i]]] <- V
     gc(verbose = FALSE)
   }
   return(object)
 }
  
}

check_one <- function(x){
  x.uni <- unique(x)
  out <- c()
  for(i in x.uni){
    idx <- which(x == i)
    if(length(idx) == 1){
      out <- c(out,i)
    }
  }
  return(out)
}

check_aga <- function(query,ref){
  query <- unique(query)
  out <- c()
  for(i in query){
    idx <- which(ref == i)
    if(length(idx) == 1){
      out <- c(out,i)
    }
  }
  return(out)
}

Run.Palette2 <- function(object,
                         nn = 2,
                         modal.norm = TRUE,
                         weight = NULL,
                         do.scaler = TRUE,
                         scaler.k = 50,
                         do.var = TRUE,
                         var.k = 20,
                         lambda = 0.5,
                         nDims = 20,
                         supervised = FALSE,
                         group = NULL,
                         origin.infer = FALSE,
                         center = T,
                         scale = T,
                         nn.k = 30,
                         k = 5,
                         nn.method = "annoy",
                         annoy.metric = "euclidean",
                         L2Norm = TRUE,
                         emb_L2 = FALSE){
  
  require(igraph)
  
  g <- dataLink(object)
  modal <- object@Data.index[["assay.idx"]]$modal
  sample <- object@Data.index[["assay.idx"]]$sample
  modal.uni <- unique(modal)
  
  # 对于每个批次，寻找它的缺失模态
  L <- length(object@Assays)
  sample.uni <- unique(sample)
  ss.idx <- c() # 记录推断数据的批次索引
  mm.idx <- c() # 记录推断数据的模态索引
  
  CN.list <- list()
  for(i in sample.uni){
    idx <- which(sample == i)[1]
    CN.list[[i]] <- colnames(object@Assays[[idx]]@data)
  }
  Tname <- Reduce(c,CN.list)
  
  weight_M_list <- list()
  real_modal_mat <- list()# 计算每个模态真实的低维表达数据
  emb.list <- list()
  for(i in modal.uni){
    idx <- which(modal == i)
    if(length(idx) > 1){
      
      if(length(idx) <= length(sample.uni)){
        d <- Reduce(c,lapply(idx,function(x){
          x <- colnames(object@Assays[[x]]@data)
          x
        }))
      }else{
        next
      }
      
    }else{
      d <- colnames(object@Assays[[idx]]@data)
    }
    
    M <- Matrix::Matrix(0, nrow = length(d), ncol = length(Tname), sparse = TRUE,
                        dimnames = list(d,Tname))
    M <- as(as(M, "generalMatrix"), "CsparseMatrix")
    M_sub <- Matrix::Matrix(0, nrow = length(d), ncol = length(d), sparse = TRUE,
                            dimnames = list(d,d))
    M_sub<- as(as(M_sub, "generalMatrix"), "CsparseMatrix")
    diag(M_sub) <- 1
    M[d,d] <- M_sub
    weight_M_list[[i]] <- M
    rm(M,M_sub,d)
    
    if(length(idx) > 1){
      real_modal_mat[[i]] <- Reduce(cbind,lapply(idx,function(x){
        x <- object@Assays[[x]]@embedding
        x
      }))
      
    }else{
      real_modal_mat[[i]] <- object@Assays[[idx]]@embedding
    }
    
    tmp_emb <- matrix(0,nrow(real_modal_mat[[i]]),length(Tname))
    rownames(tmp_emb) <- rownames(real_modal_mat[[i]])
    colnames(tmp_emb) <- Tname
    tmp_emb[,colnames(real_modal_mat[[i]])] <- real_modal_mat[[i]]
    emb.list[[i]] <- tmp_emb
    
  }
  
  for(j in sample.uni){
    idxj <- which(sample == j)
    loss.modal <- setdiff(modal.uni,modal[idxj])
    if(length(loss.modal) == 0){
      next
    }
    # remain.modal <- intersect(modal.uni,modal[idxj])
    for(jj in loss.modal){
      message(paste0("Performing inference for modal \"",jj,"\" of  sample \"",j,"\" in embedding space."))
      # 从该批次出发，到达缺失模态jj的最短路径
      paths <- lapply(all_shortest_paths(g,from = j,to = jj)[["res"]],names)
      
      # 在三个模态的情况下，缺失一个模态，必然有直接路径，缺失两个模态，至少其中一个缺失模态有直接路径
      # 如果有多个直接路径，则计算皮尔森相似系数，选择每个细胞相似系数最大的模态数据进行推断
      if(length(paths[[1]]) == 4){
        # 全为直接路径，判断是否一致
        start.m <- Reduce(c,lapply(paths, function(x) x[2]))#找到各最短路径各自的起始模态
        start.m.uni <- unique(start.m)
        if(length(start.m.uni) == 1){
          # 说明只存在一种路径，可能横跨了多个批次
          # 找到该批次出发模态数据
          start.data <- object@Assays[[intersect(which(sample == j),which(modal == start.m.uni))]]@embedding
          if(L2Norm){
            start.data <- t(L2Norm1(t(object@Assays[[intersect(which(sample == j),which(modal == start.m.uni))]]@embedding)))
          }else{
            start.data <- object@Assays[[intersect(which(sample == j),which(modal == start.m.uni))]]@embedding
          }
          # 下面判断是否横跨多个批次
          if(length(paths) > 1){
            # 路径数量大于1，说明横跨了多个批次
            
            # 找到其他批次途径模态数据
            through.data <- Reduce(cbind, lapply(1:length(paths),function(x){
              tmp.s <- paths[[x]][3]
              x <- object@Assays[[intersect(which(sample == tmp.s),which(modal == start.m.uni))]]@embedding
              x
            }))
            # 找到其他批次到达模态数据
            
            arrive.data <- Reduce(cbind, lapply(1:length(paths),function(x){
              tmp.s <- paths[[x]][3]
              x <- object@Assays[[intersect(which(sample == tmp.s),which(modal == jj))]]@embedding
              x
            }))
            
            arrive_name <- colnames(arrive.data)
          }else{
            #路径数量等于1，说明只途径一个批次
            #找到其他批次途径模态数据
            through.data <- object@Assays[[intersect(which(sample == paths[[1]][3]),which(modal == start.m.uni))]]@embedding
            # 找到其他批次到达模态数据
            arrive.data <- object@Assays[[intersect(which(sample == paths[[1]][3]),which(modal == jj))]]@embedding
            arrive_name <- colnames(arrive.data)
          }
          target <- real_modal_mat[[jj]]
          
          if(L2Norm){
            start.data <- t(L2Norm1(t(start.data)))
            through.data <- t(L2Norm1(t(through.data)))
            # arrive.data <- t(L2Norm1(t(arrive.data)))
          }
          
          # 得到出一个加权矩阵用于计算推断模态
          infer_result <- modal_infer_one_one(start.data,
                                              through.data,
                                              nn,
                                              L2Norm,
                                              arrive.data,
                                              target,
                                              do.scaler,
                                              scaler.k,
                                              do.var,
                                              var.k)
          weight_M <- infer_result[["M"]]
          colnames(weight_M) <- arrive_name
          rownames(weight_M) <- colnames(start.data)
          weight_M_list[[jj]][arrive_name,colnames(start.data)] <- t(weight_M)
          infer.data <- infer_result[[2]]
          rownames(infer.data) <- rownames(target)
          colnames(infer.data) <- colnames(start.data)
          emb.list[[jj]][,colnames(start.data)] <- infer.data

          ss.idx <- c(ss.idx,j)
          mm.idx <- c(mm.idx,jj)
          
          if(origin.infer){
            arrive_o <- Reduce(cbind, lapply(1:length(paths),function(x){
              tmp.s <- paths[[x]][3]
              x <- object@Assays[[intersect(which(sample == tmp.s),which(modal == jj))]]@data# [features,]
              x
            }))
            
            infer_o <-  arrive_o%*%t(weight_M)
            rownames(infer_o) <- rownames(arrive_o)
            colnames(infer_o) <- colnames(start.data)
            tmp1 <- Create.Musa.Assay(infer_o,
                                      filter = FALSE,
                                      min.cells = 0, 
                                      min.features = 0, 
                                      modal = jj,
                                      sample = j,
                                      slot = "data",
                                      sparce = TRUE)
            
            object <- Add.Assay(tmp1,object,check = F,class = "infer")
            rm(tmp1,infer_o,arrive_o)
          }
          rm(start.data,
             through.data,
             target,
             arrive_name,
             infer_result,
             weight_M,
             infer.data)
          gc(verbose = FALSE)
          
        }else{
          
          #说明存在不同类型的路径，途经不同模态
          start.list <- lapply(start.m.uni, function(x){
            idx1 <- which(sample == j)
            idx2 <- which(modal %in% x)
            idx <- intersect(idx1,idx2)
            x <- object@Assays[[idx]]@embedding
            x
          })
          
          mm <- Reduce(c,lapply(paths, function(x) x[2]))#得到每条路径的途径模态
          through.list <- lapply(start.m.uni, function(x){
            idx <- which(mm == x)
            if(length(idx) == 1){
              x <- object@Assays[[intersect(which(modal == x),which(sample == paths[[idx]][3]))]]@embedding
            }else{
              x <- Reduce(cbind,lapply(idx,function(y){
                y <- object@Assays[[intersect(which(modal == x),which(sample == paths[[y]][3]))]]@embedding
              }))
            }
            x
          })
          
          if(L2Norm){
            start.list <- lapply(start.list,function(x){
              x <- t(L2Norm1(t(x)))
              x
            })
            
            through.list <- lapply(through.list,function(x){
              x <- t(L2Norm1(t(x)))
              x
            })
          }
          
          arrive.data <- Reduce(cbind,lapply(start.m.uni, function(x){
            idx <- which(mm == x)
            if(length(idx) == 1){
              x <- object@Assays[[intersect(which(modal == jj),which(sample == paths[[idx]][3]))]]@embedding
            }else{
              x <- Reduce(cbind,lapply(idx,function(y){
                y <- object@Assays[[intersect(which(modal == jj),which(sample == paths[[y]][3]))]]@embedding
              }))
            }
            x
          }))
          
          arrive_name <- colnames(arrive.data)
          target <- real_modal_mat[[jj]]
          infer_result <- modal_infer_one_multi(start.list,
                                                through.list,
                                                nn,
                                                L2Norm,
                                                arrive.data,
                                                target,
                                                do.scaler,
                                                scaler.k,
                                                do.var,
                                                var.k)
          
          weight_M <- infer_result[["M"]]
          colnames(weight_M) <- arrive_name
          rownames(weight_M) <- colnames(start.list[[1]])
          weight_M_list[[jj]][arrive_name,colnames(start.list[[1]])] <- t(weight_M)
          infer.data <- infer_result[[2]]
          rownames(infer.data) <- rownames(target)
          colnames(infer.data) <- colnames(start.list[[1]])
          emb.list[[jj]][,colnames(start.list[[1]])] <- infer.data
          
          ss.idx <- c(ss.idx,j)
          mm.idx <- c(mm.idx,jj)
          
          if(origin.infer){
            
            arrive_o <- Reduce(cbind,lapply(start.m.uni, function(x){
              idx <- which(mm == x)
              if(length(idx) == 1){
                x <- object@Assays[[intersect(which(modal == jj),which(sample == paths[[idx]][3]))]]@data# [features,]
              }else{
                x <- Reduce(cbind,lapply(idx,function(y){
                  y <- object@Assays[[intersect(which(modal == jj),which(sample == paths[[y]][3]))]]@data# [features,]
                }))
              }
              x
            }))
            
            infer_o <-  arrive_o%*%t(weight_M)
            rownames(infer_o) <- rownames(arrive_o)
            colnames(infer_o) <- colnames(start.list[[1]])
            tmp1 <- Create.Musa.Assay(infer_o,
                                      filter = FALSE,
                                      min.cells = 0, 
                                      min.features = 0, 
                                      modal = jj,
                                      sample = j,
                                      slot = "data",
                                      sparce = TRUE)
            
            object <- Add.Assay(tmp1,object,check = F,class = "infer")
            rm(tmp1,infer_o,arrive_o)
          }
          rm(start.list,
             through.list,
             target,
             arrive.data,
             arrive_name,
             infer_result,
             weight_M,
             infer.data)
          gc(verbose = FALSE) 
        }
        
      }else{
        # 全为间接路径，只存在一种路径，可能横跨多个批次(在三模态马赛克成立，大于三模态则不成立)
        # 检查路径途径的模态，判断有几种间接路径（if length(paths) > 1）
        if(length(paths) == 1){
          #只有一条路径
          start.data <- object@Assays[[intersect(which(sample == j),which(modal == paths[[1]][2]))]]@embedding
          l <- length(paths[[1]])#最短路径的长度一定相等
          l.m <- seq(l/2)*2#取各路径经过的模态信息，首先创建模态索引
          
          road.list <- lapply(1:(length(l.m)-1),function(x){
            ll <- list()
            ll[[1]] <- object@Assays[[intersect(which(sample == paths[[1]][l.m[x]+1]),which(modal == paths[[1]][l.m[x]]))]]@embedding
            ll[[2]] <- object@Assays[[intersect(which(sample == paths[[1]][l.m[x]+1]),which(modal == paths[[1]][l.m[x+1]]))]]@embedding
            x <- ll
          })
          arrive.data <- object@Assays[[intersect(which(sample == paths[[1]][l.m[length(l.m)]-1]),which(modal == paths[[1]][l.m[length(l.m)]]))]]@embedding
          arrive_name <- colnames(arrive.data)
          
          if(L2Norm){
            start.data <- t(L2Norm1(t(start.data)))
            road.list <- lapply(road.list,function(x){
              x <- lapply(x,function(y){
                y <- t(L2Norm1(t(y)))
                y
              })
              x
            })
          }
          target <- real_modal_mat[[jj]]
          infer_result <- modal_infer_multi_one(start.data,
                                                road.list,
                                                nn,
                                                L2Norm,
                                                arrive.data,
                                                target,
                                                do.scaler,
                                                scaler.k,
                                                do.var,
                                                var.k)
          
          weight_M <- infer_result[["M"]]
          colnames(weight_M) <- arrive_name
          rownames(weight_M) <- colnames(start.data)
          weight_M_list[[jj]][arrive_name,colnames(start.data)] <- t(weight_M)

          infer.data <- infer_result[[2]]

          rownames(infer.data) <- rownames(target)
          colnames(infer.data) <- colnames(start.data)
          emb.list[[jj]][,colnames(start.data)] <- infer.data
          

          ss.idx <- c(ss.idx,j)
          mm.idx <- c(mm.idx,jj)
          
          if(origin.infer){
            arrive_o <- object@Assays[[intersect(which(sample == paths[[1]][l.m[length(l.m)]-1]),which(modal == paths[[1]][l.m[length(l.m)]]))]]@data# [features,]
            infer_o <- arrive_o%*%t(weight_M)
            rownames(infer_o) <- rownames(arrive_o)
            colnames(infer_o) <- colnames(start.data)
            tmp1 <- Create.Musa.Assay(infer_o,
                                      filter = FALSE,
                                      min.cells = 0, 
                                      min.features = 0, 
                                      modal = jj,
                                      sample = j,
                                      slot = "data",
                                      sparce = TRUE)
            
            object <- Add.Assay(tmp1,object,check = F,class = "infer")
            rm(tmp1,infer_o,arrive_o)
          }
          
          rm(start.data,
             road.list,
             target,
             arrive.data,
             arrive_name,
             infer_result,
             weight_M,
             infer.data)
          gc(verbose = FALSE)
        }else{
          #存在多条路径，判断是否是同一种
          l <- length(paths[[1]])#最短路径的长度一定相等
          l.m <- seq(l/2)*2#取各路径经过的模态信息，首先创建模态索引
          through.m <- lapply(paths, function(x) x[l.m])
          through.m.uni <- unique(through.m)#是否同一种路径
          
          if(length(through.m.uni) == 1){
            #是同一种路径，只需将批次合并
            
            if(L2Norm){
              start.data <- t(L2Norm1(t(object@Assays[[intersect(which(sample == j),which(modal == paths[[1]][2]))]]@embedding)))
            }else{
              start.data <- object@Assays[[intersect(which(sample == j),which(modal == paths[[1]][2]))]]@embedding
            }
            
            if(L2Norm){
              road.list <- lapply(1:(length(l.m)-1),function(x){
                ll <- list()
                ll[[1]] <- Reduce(cbind,lapply(1:length(paths),function(y){
                  y <- t(L2Norm1(t(object@Assays[[intersect(which(sample == paths[[y]][l.m[x]+1]),which(modal == paths[[y]][l.m[x]]))]]@embedding)))
                }))
                ll[[1]] <- ll[[1]][,!duplicated(colnames(ll[[1]]))]
                ll[[2]] <- Reduce(cbind,lapply(1:length(paths),function(y){
                  y <- t(L2Norm1(t(object@Assays[[intersect(which(sample == paths[[y]][l.m[x]+1]),which(modal == paths[[y]][l.m[x+1]]))]]@embedding)))
                }))
                ll[[2]] <- ll[[2]][,!duplicated(colnames(ll[[2]]))]
                ll[[1]] <- ll[[1]][,colnames(ll[[2]])]
                x <- ll
              })
            }else{
              road.list <- lapply(1:(length(l.m)-1),function(x){
                ll <- list()
                ll[[1]] <- Reduce(cbind,lapply(1:length(paths),function(y){
                  y <- object@Assays[[intersect(which(sample == paths[[y]][l.m[x]+1]),which(modal == paths[[y]][l.m[x]]))]]@embedding
                }))
                ll[[1]] <- ll[[1]][,!duplicated(colnames(ll[[1]]))]
                ll[[2]] <- Reduce(cbind,lapply(1:length(paths),function(y){
                  y <- object@Assays[[intersect(which(sample == paths[[y]][l.m[x]+1]),which(modal == paths[[y]][l.m[x+1]]))]]@embedding
                }))
                ll[[2]] <- ll[[2]][,!duplicated(colnames(ll[[2]]))]
                ll[[1]] <- ll[[1]][,colnames(ll[[2]])]
                x <- ll
              })
            }
            
            target <- real_modal_mat[[jj]]
            arrive.data <- Reduce(cbind,lapply(1:length(paths),function(y){
              y <- object@Assays[[intersect(which(sample == paths[[y]][l.m[length(l.m)]-1]),which(modal == paths[[y]][l.m[length(l.m)]]))]]@embedding
            }))
            arrive.data <- arrive.data[,!duplicated(colnames(arrive.data))]
            arrive_name <- colnames(arrive.data)
            target <- target[,arrive_name]
            infer_result <- modal_infer_multi_one(start.data,
                                                  road.list,
                                                  nn,
                                                  L2Norm,
                                                  arrive.data,
                                                  target,
                                                  do.scaler,
                                                  scaler.k,
                                                  do.var,
                                                  var.k)
            
            weight_M <- infer_result[["M"]]
            colnames(weight_M) <- arrive_name
            rownames(weight_M) <- colnames(start.data)
            weight_M_list[[jj]][arrive_name,colnames(start.data)] <- t(weight_M)
            infer.data <- infer_result[[2]]
            rownames(infer.data) <- rownames(target)
            colnames(infer.data) <- colnames(start.data)
            emb.list[[jj]][,colnames(start.data)] <- infer.data
            
            ss.idx <- c(ss.idx,j)
            mm.idx <- c(mm.idx,jj)
            
            if(origin.infer){
              arrive_o <- object@Assays[[intersect(which(sample == paths[[y]][l.m[length(l.m)]-1]),which(modal == paths[[y]][l.m[length(l.m)]]))]]@data# [features,]
              infer_o <- arrive_o%*%t(weight_M)
              rownames(infer_o) <- rownames(arrive_o)
              colnames(infer_o) <- colnames(start.data)
              tmp1 <- Create.Musa.Assay(infer_o,
                                        filter = FALSE,
                                        min.cells = 0, 
                                        min.features = 0, 
                                        modal = jj,
                                        sample = j,
                                        slot = "data",
                                        sparce = TRUE)
              
              object <- Add.Assay(tmp1,object,check = F,class = "infer")
              rm(tmp1,infer_o,arrive_o)
            }
            rm(start.data,
               road.list,
               target,
               arrive.data,
               arrive_name,
               infer_result,
               weight_M,
               infer.data)
            gc(verbose = FALSE)
            
          }else{
            # 存在多种路径，分别考虑
            # 首先将同一种路径合并
            # 得到不同路径的索引
            idx <- lapply(through.m.uni,function(x){
              tmp <- c()
              for(i in 1:length(through.m)){
                if(through.m[[i]] == x){
                  tmp <- c(tmp,i)
                }
              }
              x <- tmp
            })
            names(idx) <- through.m.uni
            # 对不同路径的批次进行合并
            
            if(L2Norm){
              start.list <- lapply(through.m.uni, function(x){
                x <- t(L2Norm1(t(object@Assays[[intersect(which(modal == x[1]),which(sample == j))]]@embedding)))
              })
            }else{
              start.list <- lapply(through.m.uni, function(x){
                x <- object@Assays[[intersect(which(modal == x[1]),which(sample == j))]]@embedding
              })
            }
            
            if(L2Norm){
              road.multi.list <- lapply(idx,function(z){
                
                z <- lapply(1:(length(l.m)-1),function(x){
                  ll <- list()
                  ll[[1]] <- Reduce(cbind,lapply(1:length(z),function(y){
                    y <- t(L2Norm1(t(object@Assays[[intersect(which(sample == paths[[z[y]]][l.m[x]+1]),which(modal == paths[[z[y]]][l.m[x]]))]]@embedding)))
                  }))
                  ll[[1]] <- ll[[1]][,!duplicated(colnames(ll[[1]]))]
                  ll[[2]] <- Reduce(cbind,lapply(1:length(z),function(y){
                    y <- t(L2Norm1(t(object@Assays[[intersect(which(sample == paths[[z[y]]][l.m[x]+1]),which(modal == paths[[z[y]]][l.m[x+1]]))]]@embedding)))
                  }))
                  ll[[2]] <- ll[[2]][,!duplicated(colnames(ll[[2]]))]
                  ll[[1]] <- ll[[1]][,colnames(ll[[2]])]
                  x <- ll
                })
                
              })
            }else{
              road.multi.list <- lapply(idx,function(z){
                
                z <- lapply(1:(length(l.m)-1),function(x){
                  ll <- list()
                  ll[[1]] <- Reduce(cbind,lapply(1:length(z),function(y){
                    y <- object@Assays[[intersect(which(sample == paths[[z[y]]][l.m[x]+1]),which(modal == paths[[z[y]]][l.m[x]]))]]@embedding
                  }))
                  ll[[1]] <- ll[[1]][,!duplicated(colnames(ll[[1]]))]
                  ll[[2]] <- Reduce(cbind,lapply(1:length(z),function(y){
                    y <- object@Assays[[intersect(which(sample == paths[[z[y]]][l.m[x]+1]),which(modal == paths[[z[y]]][l.m[x+1]]))]]@embedding
                  }))
                  ll[[2]] <- ll[[2]][,!duplicated(colnames(ll[[2]]))]
                  ll[[1]] <- ll[[1]][,colnames(ll[[2]])]
                  x <- ll
                })
                
              })
            }
            target <- real_modal_mat[[jj]]
            
            arrive.data <- Reduce(cbind,lapply(road.multi.list,function(x){
              L <- length(x)
              x <- x[[l]][[2]]
              return(x)
            }))
            arrive.data <- arrive.data[,!duplicated(colnames(arrive.data))]
            arrive_name <- colnames(arrive.data)
            target <- target[,colnames(arrive.data)]
            infer_result <- modal_infer_multi_multi(start.list,
                                                    road.multi.list,
                                                    nn,
                                                    L2Norm,
                                                    arrive.data,
                                                    target,
                                                    do.scaler,
                                                    scaler.k,
                                                    do.var,
                                                    var.k)
            
            weight_M <- infer_result[["M"]]
            colnames(weight_M) <- arrive_name
            if(length(unique(arrive_name)) < length(arrive_name)){
              unique_arrive_name <- unique(arrive_name)
              # 创建映射矩阵
              map_matrix <- sapply(unique_arrive_name, function(x) arrive_name == x)
              map_matrix <- matrix(as.numeric(map_matrix), nrow = length(arrive_name))
              weight_M <- weight_M %*% map_matrix
              colnames(weight_M) <- unique_arrive_name
            }
            
            rownames(weight_M) <- colnames(start.list[[1]])
            weight_M_list[[jj]][arrive_name,colnames(start.list[[1]])] <- t(weight_M)

            infer.data <- infer_result[[2]]

            colnames(infer.data) <- colnames(start.list[[1]])
            emb.list[[jj]][,colnames(start.list[[1]])] <- infer.data

            ss.idx <- c(ss.idx,j)
            mm.idx <- c(mm.idx,jj)
            if(origin.infer){
              # features <- object@Co.Features[[jj]]
              arrive_o <- Reduce(cbind,lapply(idx,function(x){
                x <- Reduce(cbind,lapply(x,function(y){
                  y <- object@Assays[[intersect(which(sample == paths[[y]][max(l.m)-1]),which(modal == jj))]]@data# [features,]
                  return(y)
                }))
                return(x)
              }))
              infer_o <- arrive_o%*%t(weight_M)
              rownames(infer_o) <- rownames(arrive_o)
              colnames(infer_o) <- colnames(start.list[[1]])
              tmp1 <- Create.Musa.Assay(infer_o,
                                        filter = FALSE,
                                        min.cells = 0, 
                                        min.features = 0, 
                                        modal = jj,
                                        sample = j,
                                        slot = "data",
                                        sparce = TRUE)
              
              object <- Add.Assay(tmp1,object,check = F,class = "infer")
              rm(tmp1,infer_o,arrive_o)
            }
            rm(start.list,
               road.multi.list,
               target,
               arrive.data,
               arrive_name,
               infer_result,
               weight_M,
               infer.data)
            gc(verbose = FALSE)
            
          }
        }
      }
    }
  }

  emb.list <- lapply(emb.list, function(x){
    x <- ScaleData1(object = x,
                    features = rownames(x),
                    do.scale = scale,
                    do.center = center,
                    scale.max = 1000000,
                    block.size = 1000,
                    min.cells.to.block = 3000,
                    verbose = FALSE)
    x
  })
  
  if(modal.norm){
    C.list <- lapply(emb.list,function(x){
      x <- sum(sqrt(diag(tcrossprod(x))))
      x
    })
    maxc <- max(unlist(C.list))
    emb.list <- lapply(1:length(emb.list),function(x){
      x <- emb.list[[x]]*maxc/C.list[[x]]
      x
    })
  }
  
  names(emb.list) <- modal.uni
  for(i in modal.uni){
    rownames(emb.list[[i]]) <- paste0(i,rownames(emb.list[[i]]))
  }
  
  if(!is.null(weight)){
    embname <- names(emb.list)
    if(any(unlist(weight) == 0)){
      weight <- weight[names(emb.list)]
      idx <- which(unlist(weight) == 0)
      emb.list <- emb.list[-idx]
      weight <- weight[-idx]
      embname <- names(emb.list)
    }
    emb.list <- lapply(modal.uni,function(x){
      x <- emb.list[[x]]*weight[[x]]
      x
    })
    
    names(emb.list) <- embname
  }
  
  
  emb <- Reduce(rbind,emb.list)
  
  # 以前的代码，使用这个删除if(supervised){}elses{这里的部分}中else的部分，恢复if(supervised){}else{}下面的部分
  # SNN <- as(FindNeighbors1(t(emb),
  #                          k.param = nn.k,
  #                          nn.method = nn.method,
  #                          prune.SNN = 0,
  #                          annoy.metric = annoy.metric,
  #                          verbose = FALSE)$snn,"dgCMatrix")
  # SNN <- SNN[Tname,Tname]
  # emb <- emb[,Tname]
  # index.list <- list()
  # for(i in 1:length(CN.list)){
  #   index.list[[i]] <- match(CN.list[[i]],Tname)
  # }
  if(supervised && is.null(group)){
    stop('Must selcet meta data under supervised case.')
  }
  if(supervised){
    if(length(intersect(group,colnames(object@Meta))) == 0){
      stop('Group information does not exist.')
    }
    meta_cluster <- as.character(object@Meta[,group])
    
    ####### 缺失细胞类型
    idx <- which(!is.na(meta_cluster))
    names(meta_cluster) <- object@Meta[,'cell.name'][idx]
    emb_sub <- emb[,object@Meta[,'cell.name'][idx]]
    SNN <- as(FindNeighbors1(t(emb_sub),
                             k.param = nn.k,
                             nn.method = nn.method,
                             prune.SNN = 0,
                             annoy.metric = annoy.metric,
                             verbose = FALSE)$snn,"dgCMatrix")
    SNN <- SNN[object@Meta[,'cell.name'][idx],object@Meta[,'cell.name'][idx]]
    
    index.list <- list()
    meta_list <- list()
    for(i in 1:length(CN.list)){
      tmp <- intersect(CN.list[[i]],object@Meta[,'cell.name'][idx])
      index.list[[i]] <- match(tmp,object@Meta[,'cell.name'][idx])
      meta_list[[i]] <- meta_cluster[index.list[[i]]]
    }
    SNN <- filter_SNN(SNN,index.list,meta_list)
    colnames(SNN) <- rownames(SNN) <- object@Meta[,'cell.name'][idx]
    
    
    K_list <- K_SNN(SNN,index.list,k,lambda)
    K <- K_list[["K"]]
    rownames(K) <- colnames(K) <- object@Meta[,'cell.name'][idx]
    B <- Matrix::Matrix(0,nrow = nrow(SNN),ncol = ncol(SNN),dimnames = list(object@Meta[,'cell.name'][idx],object@Meta[,'cell.name'][idx]))
    diag(B) <- K_list[["B"]]
    Zname <- which(rowSums(K) == 0)
    diag(K)[Zname] <- K_list[["value"]]
    
    mat <- Getmat(emb_sub,K,B)
    rownames(mat) <- colnames(mat) <- rownames(emb_sub)
    
    Subspace <- suppressWarnings(RSpectra::eigs_sym(mat,k = nDims,which = "LA"))[["vectors"]] %>% as.matrix()
    
    #######
    
    # 以前的代码
    # names(meta_cluster) <- object@Meta[,'cell.name']
    # meta_cluster <- as.character(meta_cluster[Tname])
    # meta_list <- list()
    # for(i in 1:length(CN.list)){
    #   meta_list[[i]] <- meta_cluster[index.list[[i]]]
    # }
    # 
    # SNN <- filter_SNN(SNN,index.list,meta_list)
    # colnames(SNN) <- rownames(SNN) <- Tname
  }else{
    SNN <- as(FindNeighbors1(t(emb),
                             k.param = nn.k,
                             nn.method = nn.method,
                             prune.SNN = 0,
                             annoy.metric = annoy.metric,
                             verbose = FALSE)$snn,"dgCMatrix")
    SNN <- SNN[Tname,Tname]
    emb <- emb[,Tname]
    index.list <- list()
    for(i in 1:length(CN.list)){
      index.list[[i]] <- match(CN.list[[i]],Tname)
    }
    
    K_list <- K_SNN(SNN,index.list,k,lambda)
    K <- K_list[["K"]]
    rownames(K) <- colnames(K) <- Tname
    B <- Matrix::Matrix(0,nrow = nrow(SNN),ncol = ncol(SNN),dimnames = list(Tname,Tname))
    diag(B) <- K_list[["B"]]
    Zname <- which(rowSums(K) == 0)
    diag(K)[Zname] <- K_list[["value"]]

    mat <- Getmat(emb,K,B)
    rownames(mat) <- colnames(mat) <- rownames(emb)

    Subspace <- suppressWarnings(RSpectra::eigs_sym(mat,k = nDims,which = "LA"))[["vectors"]] %>% as.matrix()
  }
  
  # 以前的代码，使用这个删除if(supervised){}elses{这里的部分}中else的部分，恢复if(supervised)上面的部分
  # K_list <- K_SNN(SNN,index.list,k,lambda)
  # K <- K_list[["K"]]
  # rownames(K) <- colnames(K) <- Tname
  # B <- Matrix::Matrix(0,nrow = nrow(SNN),ncol = ncol(SNN),dimnames = list(Tname,Tname))
  # diag(B) <- K_list[["B"]]
  # Zname <- which(rowSums(K) == 0)
  # diag(K)[Zname] <- K_list[["value"]]
  # 
  # mat <- Getmat(emb,K,B)
  # rownames(mat) <- colnames(mat) <- rownames(emb)
  # 
  # Subspace <- suppressWarnings(RSpectra::eigs_sym(mat,k = nDims,which = "LA"))[["vectors"]] %>% as.matrix()
  
  embedding <- crossprod(as.matrix(Reduce(rbind,emb.list)), Subspace)
  if(emb_L2){
    embedding <- L2Norm1(embedding)
  }
  rownames(embedding) <- Tname
  colnames(embedding) <- colnames(Subspace) <- paste0("dim_",seq(ncol(embedding)))
  object@Int.result[["bind"]] <- embedding
  
  return(object)
}

Trans.Knowledge <- function(object,
                            query,
                            ref,
                            ref.label = NULL,
                            k = 5){
  
  
  # object:计算对象
  # query：查询批次
  # ref:参考批次
  # ref.label:一个一列的数据框，行名为细胞名，列中内容为待转移的标签
  # k:使用多少k近邻对标签进行预测
  
  data <- object@Int.result[["bind"]]
  
  
  
  
  pred_labels_signal = class::knn(t(res[, which(Meta$dataset == "reference")]),
                                  t(res[, which(Meta$dataset == "query")]),
                                  as.character(Meta$SubClass[which(Meta$dataset == "reference")]),
                                  k = 5, prob = TRUE) %>% as.character()
  
}

############## reference-based integration code ##############
library(reticulate)
library(batchelor)
library(harmony)
library(Seurat)
library(Signac)
library(scater)
library(irlba)
library(BiocNeighbors)
library(SingleCellExperiment)
library(Matrix)
library(umap)
library(dplyr)
library(Rcpp)

get_representation <- function(query,
                               ref,
                               latent,
                               Bi_sPCA = TRUE,
                               dims = 30,
                               lambda = 0.5,
                               do.scale = TRUE,
                               do.center = TRUE,
                               nn.k = 20,# knn snn参数
                               nn.method = "annoy",# knn snn参数
                               annoy.metric = "euclidean",
                               prune.SNN = 1/15,
                               verbose = TRUE){
  
  # ref和query经过批次整合后输入函数，得到query在马赛克整合嵌入上的表示
  # 本函数输出结果再进行一步整合即完成ref-query 映射
  
  # query: query data matrix after batch integrate,feature by cell
  # ref: reference data matrix after batch integrate,feature by cell
  # latent: reference的已完成整合的低维嵌入,feature by cell
  # Bi_sPCA：是否执行Bi-sPCA
  # lambda: 参数
  
  
  if(all(colnames(ref) %in% colnames(latent))){
    latent <- latent[,colnames(ref)]
  }else{
    stop('Unmatched reference and latent')
  }
  
  if(nrow(query) != nrow(ref)){
    stop('Unmatched dim between reference and query')
  }
  
  if(Bi_sPCA){
    if(dims >= nrow(ref)){
      Bi_sPCA = FALSE
      message('Input dimension is lower than the preset reduced dimension, set Bi_sPCA to FALSE.')
    }
    
    if(lambda < 0){
      lambda = 0
      message('lambda should be between 0 and 1, set lambda to 0.')
    }
    
    if(lambda > 1){
      lambda = 1
      message('lambda should be between 0 and 1, set lambda to 1.')
    }
  }
  
  latent <- latent[,colnames(ref)]
  
  if(Bi_sPCA){
    K1 <- FindNeighbors1(t(as.matrix(latent)),
                         k.param = nn.k,
                         nn.method = nn.method,
                         prune.SNN = prune.SNN,
                         annoy.metric = annoy.metric,
                         verbose = verbose)$snn
    
    degrees = colSums(K1)
    degrees[which(degrees == 0)] <- 1e-8
    D <- Matrix::Diagonal(x = 1 / sqrt(degrees))
    K1 <- D %*% K1 %*% D
    
    K2 <- FindNeighbors1(t(as.matrix(ref)),
                         k.param = nn.k,
                         nn.method = nn.method,
                         prune.SNN = prune.SNN,
                         annoy.metric = annoy.metric,
                         verbose = verbose)$snn
    
    degrees = colSums(K2)
    degrees[which(degrees == 0)] <- 1e-8
    D <- Matrix::Diagonal(x = 1 / sqrt(degrees))
    K2 <- D %*% K2 %*% D
    
    rm(degrees,D)
    
    obj <- cbind(query,ref)
    rownames(obj) <- paste0('dim',1:nrow(obj))
    rownames(ref) <- rownames(query) <- rownames(obj)  
    
    data <- ScaleData1(object = obj,
                       features = rownames(ref),
                       do.scale = do.scale,
                       do.center = do.center,
                       scale.max = 10,
                       block.size = 1000,
                       min.cells.to.block = 3000,
                       verbose = FALSE)
    rm(obj)
    
    colnames(data) <- c(colnames(query),colnames(ref))
    s_ref <- data[,colnames(ref)]
    s_query <- data[,colnames(query)]
    
    mat <- s_ref%*%(K1-lambda*K2)%*%t(s_ref)
    
    rownames(mat) <- colnames(mat) <- rownames(ref)
    Subspace <- suppressWarnings(RSpectra::eigs_sym(mat,k = dims,which = "LA"))[["vectors"]] %>% as.matrix()
    
    embedding <- t(Subspace) %*% data
    colnames(embedding) <- c(colnames(query),colnames(ref))
    rownames(embedding) <- paste0('dim',1:nrow(embedding))
    ref <- embedding[,colnames(ref)]
    query <- embedding[,colnames(query)]
    
  }else{
    ref <- as.matrix(ref)
    query <- as.matrix(query)
  }
  
  latent <- as.matrix(latent)
  
  out <- latent %*% MASS::ginv(ref) %*% query
  colnames(out) <- colnames(query)
  rownames(out) <- rownames(latent)
  return(out)
}

run_fastMNN = function(batches, 
                       cell_name, 
                       is.normalize = TRUE,
                       method = 'LogNormalize',
                       nfeatures = NULL,
                       vargs = NULL,
                       k = 20, 
                       out.npcs = NA){
  # preprocessing
  
  seurat.list <- sapply(batches, function(x){
    CreateSeuratObject(x)
  }, simplify = FALSE)
  # modified rownames
  seurat.list <- lapply(X = seurat.list, FUN = function(x) {
    # VariableFeatures(x) = rownames(batches[[1]])
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = Inf)
  })
  # VariableFeatures(batch_seurat) = vargs
  vargenes <- seurat.list[[1]]@assays[["RNA"]]@var.features
  
  if (is.normalize == TRUE) {
    if(method == 'LogNormalize'){
      seurat.list <- lapply(X = seurat.list, FUN = function(x) {
        x <- NormalizeData(x)
        x
      })
    }
    if(method == 'CLR'){
      seurat.list <- lapply(X = seurat.list, FUN = function(x) {
        x <- NormalizeData(x, normalization.method = 'CLR', margin = 2)
        x
      })
    }
    
  }
  
  if(method == 'TF-IDF'){
    obj <- CreateSeuratObject(Reduce(cbind,batches))
    obj <-  RunTFIDF(obj)
    obj[['lsi']] <- RunSVD(obj@assays[[1]]@data)
    
    data_list = lapply(1:length(batches), function(i) t(obj[['lsi']]@cell.embeddings[colnames(batches[[i]]),2:(out.npcs+1)]))
    rm(obj)
    
  }else{
    
    data_list = lapply(1:length(seurat.list), function(i) seurat.list[[i]]@assays[["RNA"]]@data[vargenes,])
    
  }
  
  rm(seurat.list)
  ##########################################################
  # run fastMNN
  t1 = Sys.time()
  out_mnn_total = do.call(batchelor::fastMNN, c(data_list, k = k, d = unname(out.npcs)))
  t2 = Sys.time()
  print(t2-t1)
  
  fastmnn_res = t(out_mnn_total@assays@data@listData[["reconstructed"]]@seed@components)
  colnames(fastmnn_res) = cell_name
  rownames(fastmnn_res) <- paste0('dim',1:nrow(fastmnn_res))
  return(fastmnn_res)
}

run_Seurat = function(batches, 
                      cell_name, 
                      is.normalize = TRUE,
                      method = 'LogNormalize',
                      nfeatures = NULL,
                      vargs = NULL,
                      k = 20, 
                      dims = NULL,
                      out.npcs = 30){
  # preprocessing
  
  seurat.list <- sapply(batches, function(x){
    CreateSeuratObject(x)
  }, simplify = FALSE)
  # modified rownames
  seurat.list <- lapply(X = seurat.list, FUN = function(x) {
    # VariableFeatures(x) = rownames(batches[[1]])
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = Inf)
  })
  # VariableFeatures(batch_seurat) = vargs
  vargenes <- seurat.list[[1]]@assays[["RNA"]]@var.features
  
  if (is.normalize == TRUE) {
    if(method == 'LogNormalize'){
      seurat.list <- lapply(X = seurat.list, FUN = function(x) {
        x <- NormalizeData(x)
      })
    }
    if(method == 'CLR'){
      seurat.list <- lapply(X = seurat.list, FUN = function(x) {
        x <- NormalizeData(x, normalization.method = 'CLR', margin = 2)
      })
    }
    
    if(method == 'TF-IDF'){
      seurat.list <- lapply(X = seurat.list, FUN = function(x) {
        x <- RunTFIDF(x)
      })
    }
  }
  
  if(method == 'TF-IDF'){
    seurat.list <- lapply(X = seurat.list, FUN = function(x) {
      # x <- RunSVD(x,reduction.name = 'lsi')
      x[['lsi']] <-  RunSVD(x@assays[[1]]@data)
      x
    })
    cell_anchors = FindIntegrationAnchors(object.list = seurat.list, 
                                          k.filter = 200,
                                          reduction = 'rlsi', 
                                          anchor.features = vargenes,
                                          dims = 2:(out.npcs+1))
    rm(seurat.list)
    comb <- CreateSeuratObject(do.call(cbind, batches))
    VariableFeatures(comb) = rownames(batches[[1]])
    comb <- RunTFIDF(comb)
    comb <- RunSVD(comb,reduction.name = 'lsi',n = out.npcs)
    gc()
    batch_correct = IntegrateEmbeddings(anchorset = cell_anchors,
                                        reductions = comb[["lsi"]],
                                        new.reduction.name = "integrated_lsi",
                                        dims.to.integrate = 1:out.npcs)
    rm(comb,cell_anchors)
    seurat_res <- t(as.matrix(batch_correct[["integrated_lsi"]]@cell.embeddings))
    
    rm(batch_correct)
  }else{
    cell_anchors = FindIntegrationAnchors(object.list = seurat.list,
                                          k.filter = 200,
                                          anchor.features = vargenes,
                                          reduction = 'cca')
    rm(seurat.list)
    if(is.na(out.npcs)){
      
      if(is.null(dims)){
        dims = 1:length(vargenes)
      }
      
      batch_correct = IntegrateData(anchorset = cell_anchors,dims = dims)
      seurat_res <- as.matrix(batch_correct@assays[["integrated"]]@data)
    }else{
      batch_correct = IntegrateData(anchorset = cell_anchors)
      DefaultAssay(batch_correct) = "integrated"
      batch_correct = ScaleData(object = batch_correct)
      batch_correct = RunPCA(object = batch_correct, npcs = out.npcs, verbose = FALSE)
      seurat_res = t(as.matrix(batch_correct@reductions$pca@cell.embeddings))
    }
  }
  
  colnames(seurat_res) = cell_name
  return(seurat_res)
}

run_Harmony = function(batches, cell_name,
                       is.normalize = TRUE,
                       method = 'LogNormalize',
                       theta = 2, 
                       out.npcs = 30){
  # preprocessing
  batch_seurat = CreateSeuratObject(do.call(cbind, batches))
  batch_info <- unlist(lapply(1:length(batches),FUN = function(x){
    x <- rep(paste0("bath",x),ncol(batches[[x]]))
    x
  }))
  batch_seurat[["batch"]] <- batch_info
  # modified rownames
  batch_seurat <- FindVariableFeatures(batch_seurat, selection.method = "vst", nfeatures = Inf)
  # VariableFeatures(batch_seurat) = rownames(batches[[1]])
  # VariableFeatures(batch_seurat) = vargs
  vargenes <- batch_seurat@assays[["RNA"]]@var.features
  
  if (is.normalize == TRUE){
    if(method == 'LogNormalize'){
      batch_seurat <- NormalizeData(batch_seurat)
      if(is.na(out.npcs)){
        # batch_seurat@assays[['RNA']]@data <- batch_seurat@assays[['RNA']]@counts
        batch_seurat@assays[['RNA']]@scale.data <- as.matrix(batch_seurat@assays[['RNA']]@counts)
        batch_seurat[['pca']] <- CreateDimReducObject(embeddings = as.matrix(t(batch_seurat@assays[['RNA']]@data)),key = 'PC_')
      }else{
        batch_seurat <- ScaleData(batch_seurat)
        batch_seurat <- RunPCA(batch_seurat,npcs = out.npcs)
      }
    }
    if(method == 'CLR'){
      batch_seurat <- NormalizeData(batch_seurat, normalization.method = 'CLR', margin = 2)
      if(is.na(out.npcs)){
        # batch_seurat@assays[['RNA']]@data <- batch_seurat@assays[['RNA']]@counts
        batch_seurat@assays[['RNA']]@scale.data <- as.matrix(batch_seurat@assays[['RNA']]@counts)
        batch_seurat[['pca']] <- CreateDimReducObject(embeddings = as.matrix(t(batch_seurat@assays[['RNA']]@data)),key = 'PC_')
      }else{
        batch_seurat <- ScaleData(batch_seurat)
        batch_seurat <- RunPCA(batch_seurat,npcs = out.npcs)
      }
    }
    
    if(method == 'TF-IDF'){
      batch_seurat <- RunTFIDF(batch_seurat)
      tmp <- RunSVD(batch_seurat@assays[[1]]@data, n = (out.npcs+1))
      batch_seurat[['pca']] <- CreateDimReducObject(embeddings = as.matrix(tmp@cell.embeddings[,2:(out.npcs+1)]),key = 'PC_')
      rm(tmp)
      # batch_seurat <- RunSVD(batch_seurat, n = (out.npcs+1),
      #                        reduction.name = 'pca')
    }
  }else{
    if(is.na(out.npcs)){
      # batch_seurat@assays[['RNA']]@data <- batch_seurat@assays[['RNA']]@counts
      batch_seurat@assays[['RNA']]@scale.data <- as.matrix(batch_seurat@assays[['RNA']]@counts)
      batch_seurat[['pca']] <- CreateDimReducObject(embeddings = as.matrix(t(batch_seurat@assays[['RNA']]@data)),key = 'PC_')
    }else{
      batch_seurat <- ScaleData(batch_seurat)
      batch_seurat <- RunPCA(batch_seurat,npcs = out.npcs)
    }
  }
  
  #run Harmony
  t1 = Sys.time()
  if(method == 'TF-IDF'){
    batch_seurat = RunHarmony(object = batch_seurat, 
                              group.by.vars = 'batch', 
                              theta = theta,
                              dims.use = 1:out.npcs,
                              # dims.use = 2:(out.npcs+1),
                              project.dim = F,
                              plot_convergence = TRUE, 
                              nclust = 50, 
                              max.iter.cluster = 100) 
  }else{
    batch_seurat = RunHarmony(object = batch_seurat, 
                              group.by.vars = 'batch', 
                              theta = theta, 
                              plot_convergence = TRUE, 
                              nclust = 50, 
                              max.iter.cluster = 100)
  }
  t2 = Sys.time()
  print(t2-t1)
  
  harmony_res = t(as.matrix(batch_seurat@reductions[["harmony"]]@cell.embeddings))
  colnames(harmony_res) = cell_name
  return(harmony_res)
}

get_ref_idx <- function(q_modal, ref_batch, ref_modal) {
  # 确保 q_modal 唯一
  q_modal <- unique(q_modal)
  
  # 存储所有匹配的批次
  matched_batches <- vector("list", length(q_modal))
  
  for (i in seq_along(q_modal)) {
    # 找到 ref_modal 中对应的索引
    tmp_idx <- which(ref_modal == q_modal[i])
    if (length(tmp_idx) == 0) {
      # 如果没有匹配值，直接返回 NULL
      return(NULL)
    }
    # 将匹配到的批次存储
    matched_batches[[i]] <- ref_batch[tmp_idx]
  }
  
  # 取所有匹配批次的并集
  # final_batches <- Reduce(union, matched_batches)
  final_batches <- Reduce(intersect, matched_batches)
  
  # 找到最终匹配批次在 ref_batch 中的索引
  ref_idx <- which(ref_batch %in% final_batches)
  
  return(ref_idx)
}

split_block <- function(modal,batch){
  uni_batch <- unique(batch)
  modal_list <- list()
  for(i in uni_batch){
    tmp_idx <- which(batch == i)
    modal_list[[i]] <- modal[tmp_idx]
  }
  
  vector_signatures <- sapply(modal_list, function(vec) paste(sort(vec), collapse = ","))
  
  unique_signatures <- unique(vector_signatures)
  out <- lapply(unique_signatures, function(sig) which(vector_signatures == sig))
  
  return(out)
}

Palette_map <- function(datalist,
                        modal,
                        batch,
                        q_or_batch,
                        modal_order,
                        method = c('fastMNN','Seurat','Harmony'),
                        is.normalize = TRUE,
                        normalize_method = c('LogNormalize','CLR','TF-IDF'),
                        method_out.npcs = 30,
                        latent,
                        Bi_sPCA = TRUE,
                        dims = 30,
                        lambda = 0.5,
                        do.scale = TRUE,
                        do.center = TRUE,
                        nn.k = 20,
                        nn.method = "annoy",
                        annoy.metric = "euclidean",
                        prune.SNN = 1/15,
                        verbose = TRUE){
  
  idx_q <- which(q_or_batch == 'query')
  q_modal <- modal[idx_q]
  q_batch <- batch[idx_q]
  query <- datalist[idx_q]
  
  idx_r <- which(q_or_batch == 'ref')
  r_modal <- modal[idx_r]
  r_batch <- batch[idx_r]
  ref <- datalist[idx_r]
  
  rm(datalist,modal,batch,q_or_batch)
  
  if(!is.null(modal_order)){
    diff <- setdiff(c(unique(q_modal),modal_order),modal_order)
    if(length(diff) > 0 ){
      modal_order <- unique(q_modal)
      message('The modal_order parameter is incorrect and has been reset')
    }
  }else{
    modal_order <- unique(q_modal)
  }
  
  if(length(method) != length(modal_order)){
    method <- rep(method[1],length(modal_order))
  }
  names(method) <- modal_order
  if(length(is.normalize) != length(modal_order)){
    is.normalize <- rep(is.normalize[1],length(modal_order))
  }
  names(is.normalize) <- modal_order
  if(length(method_out.npcs) != length(modal_order)){
    method_out.npcs <- rep(method_out.npcs[1],length(modal_order))
  }
  names(method_out.npcs) <- modal_order
  if(length(normalize_method) != length(modal_order)){
    normalize_method <- rep(normalize_method[1],length(modal_order))
  }
  names(normalize_method) <- modal_order
  
  batch_uni <- unique(q_batch)
  
  batch_block <- split_block(q_modal,q_batch)
  rep_list <- list()
  
  for(i in 1:length(batch_block)){
    
    tmp_block <- batch_uni[batch_block[[i]]]
    
    
    tmp_idx <- which(q_batch %in% tmp_block)
    tmp_modal <- q_modal[tmp_idx]
    tmp_data <- query[tmp_idx]
    
    ref_idx <- get_ref_idx(tmp_modal,r_batch,r_modal)
    if(is.null(ref_idx)){
      stop(paste0('No reference could be found that matched batch ',i,'.'))
    }
    
    tmp_ref <- ref[ref_idx]
    ref_modal <- r_modal[ref_idx]
    ref_batch <- r_batch[ref_idx]
    
    int_ref <- list()
    int_query <- list()
    
    for(j in unique(tmp_modal)){
      method_j <- method[j]
      n_method_j <- normalize_method[j]
      is.normalize_j <- is.normalize[j]
      method_out.npcs_j <- method_out.npcs[j]
      
      idx_q_j <- which(tmp_modal == j)
      idx_r_j <- which(ref_modal == j)
      data_j <- c(tmp_data[idx_q_j],tmp_ref[idx_r_j])
      cell_name_j <- unlist(sapply(data_j,colnames))
      if(method_j == "Seurat"){
        res <- run_Seurat(data_j,
                          cell_name_j,
                          is.normalize = is.normalize_j,
                          method = n_method_j,
                          out.npcs = method_out.npcs_j)
      }
      
      if(method_j == "fastMNN"){
        res <- run_fastMNN(data_j,
                           cell_name_j,
                           is.normalize = is.normalize_j,
                           method = n_method_j,
                           out.npcs = unname(method_out.npcs_j))
      }
      
      if(method_j == "Harmony"){
        res <- run_Harmony(data_j,
                           cell_name_j,
                           is.normalize = is.normalize_j,
                           method = n_method_j,
                           out.npcs = unname(method_out.npcs_j))
      }
      
      int_ref[[j]] <- res[,unlist(lapply(tmp_ref[idx_r_j],colnames))]
      int_query[[j]] <- res[,unlist(lapply(tmp_data[idx_q_j],colnames))]
    }
    
    rm(res,data_j)
    if(length(int_ref)>1){
      int_ref_mat <- Reduce(rbind,int_ref)
      int_query_mat <- Reduce(rbind,int_query)
    }else{
      int_ref_mat <- int_ref[[1]]
      int_query_mat <- int_query[[1]]
    }
    
    rm(int_ref,int_query)
    rep_list[[i]] <- get_representation(int_query_mat,
                                        int_ref_mat,
                                        latent = latent,
                                        Bi_sPCA = Bi_sPCA,
                                        dims = dims,
                                        lambda = lambda,
                                        do.scale = do.scale,
                                        do.center = do.center,
                                        nn.k = nn.k,
                                        nn.method = nn.method,
                                        annoy.metric = annoy.metric,
                                        prune.SNN = prune.SNN,
                                        verbose = verbose)
    
  }
  
  return(rep_list)
}


