## code from Liu et al. https://github.com/PYangLab/scMultiBench/tree/main/evaluation_pipelines/imputation
library(tools)
library(glue)
library(data.table)
library(limma)
library(edgeR)
library(matrixStats)

compute_smse <- function(impute_data, gt_data) {
  mse <- colMeans((impute_data - gt_data)^2)
  sd_gt <- apply(gt_data, 2, sd)
  smse <- mse / (sd_gt^2)
  return(smse)
}

doLimma_singlecore <- function(exprsMat, cellTypes, exprs_pct = 0.05){
  cellTypes <- droplevels(as.factor(cellTypes))
  tt <- list()
  for (i in 1:nlevels(cellTypes)) {
    tmp_celltype <- (ifelse(cellTypes == levels(cellTypes)[i], 1, 0))
    design <- stats::model.matrix(~tmp_celltype)
    meanExprs <- do.call(cbind, lapply(c(0,1), function(i){
      Matrix::rowMeans(exprsMat[, tmp_celltype == i, drop = FALSE])
    }))
    meanPct <- do.call(cbind, lapply(c(0,1), function(i){
      Matrix::rowSums(exprsMat[, tmp_celltype == i, drop = FALSE] > 0)/sum(tmp_celltype == i)
    }))
    #keep <- meanPct[,2] > exprs_pct
    y <- methods::new("EList")
    #y$E <- exprsMat[keep, ]
    y$E <- exprsMat
    fit <- limma::lmFit(y, design = design)
    fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
    tt[[i]] <- limma::topTable(fit, n = Inf, adjust.method = "BH", coef = 2)
    if (!is.null(tt[[i]]$ID)) {
      tt[[i]] <- tt[[i]][!duplicated(tt[[i]]$ID),]
      rownames(tt[[i]]) <- tt[[i]]$ID
    }
    tt[[i]]$meanExprs.1 <- meanExprs[rownames(tt[[i]]), 1]
    tt[[i]]$meanExprs.2 <- meanExprs[rownames(tt[[i]]), 2]
    tt[[i]]$meanPct.1 <- meanPct[rownames(tt[[i]]), 1]
    tt[[i]]$meanPct.2 <- meanPct[rownames(tt[[i]]), 2]
  }
  names(tt) <- levels(cellTypes)
  return(tt)
}

doLimma <- function(exprsMat, cellTypes, exprs_pct = 0.05){
  message("Limma multicore")
  cellTypes <- droplevels(as.factor(cellTypes))
  
  tt <- mclapply(c(1:nlevels(cellTypes)), mc.cores =1, function(i){
    tmp_celltype <- (ifelse(cellTypes == levels(cellTypes)[i], 1, 0))
    design <- stats::model.matrix(~tmp_celltype)
    meanExprs <- do.call(cbind, lapply(c(0,1), function(i){
      Matrix::rowMeans(exprsMat[, tmp_celltype == i, drop = FALSE])
    }))
    meanPct <- do.call(cbind, lapply(c(0,1), function(i){
      Matrix::rowSums(exprsMat[, tmp_celltype == i, drop = FALSE] > 0)/sum(tmp_celltype == i)
    }))
    #keep <- meanPct[,2] > exprs_pct
    y <- methods::new("EList")
    #y$E <- exprsMat[keep, ]
    y$E <- exprsMat
    fit <- limma::lmFit(y, design = design)
    fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
    temp <- limma::topTable(fit, n = Inf, adjust.method = "BH", coef = 2)
    if (!is.null(temp$ID)) {
      temp <- temp[!duplicated(temp$ID),]
      rownames(temp) <- temp$ID
    }
    temp$meanExprs.1 <- meanExprs[rownames(temp), 1]
    temp$meanExprs.2 <- meanExprs[rownames(temp), 2]
    temp$meanPct.1 <- meanPct[rownames(temp), 1]
    temp$meanPct.2 <- meanPct[rownames(temp), 2]
    return(temp)
  })
  names(tt) <- levels(cellTypes)
  return(tt)
}

doTest <- function(exprsMat, cellTypes) {
  # input must be normalised, log-transformed data
  message("t-test multicore")
  cty <- droplevels(as.factor(cellTypes))
  names(cty) <- colnames(exprsMat)
  
  tt <- mclapply(c(1:nlevels(cty)), mc.cores =3, function(i){
    
    tmp_celltype <- (ifelse(cty == levels(cty)[i], 1, 0))
    
    temp <- t(apply(exprsMat, 1, function(x) {
      x1 <- x[tmp_celltype == 0]
      x2 <- x[tmp_celltype == 1]
      
      res <- stats::t.test(x2, y=x1)
      return(c(stats=res$statistic,
               pvalue=res$p.value))
    }))
    temp <- as.data.frame(temp)
    temp$adj.pvalue <- stats::p.adjust(temp$pvalue, method = "BH")
    return(temp)
  })
  names(tt) <- levels(cty)
  return(tt)
}

make_sce <- function(expr, label){
  sce <- SingleCellExperiment(list(logcounts=expr))
  sce$celltype <- as.factor(label)
  return(sce)
}

drop_zero_or_const_rows <- function(ref_mat, test_mat) {
  # 统一为普通矩阵以便 rowSums/rowVars；保留稀疏也行，但 rowVars 需要 dense
  ref <- as.matrix(ref_mat)
  tst <- as.matrix(test_mat)
  
  # 先按“全部为 0”过滤（任一矩阵中出现全 0 都会导致相关为 NA）
  keep_nonallzero <- (rowSums(ref) != 0) & (rowSums(tst) != 0)
  
  if (!any(keep_nonallzero)) {
    return(list(ref = ref[FALSE, , drop = FALSE],
                test = tst[FALSE, , drop = FALSE]))
  }
  
  ref2 <- ref[keep_nonallzero, , drop = FALSE]
  tst2 <- tst[keep_nonallzero, , drop = FALSE]
  
  # 再过滤“方差为 0”的常数行（全 1、全常数等）
  ref_var_ok <- matrixStats::rowVars(ref2) > 0
  tst_var_ok <- matrixStats::rowVars(tst2) > 0
  keep <- ref_var_ok & tst_var_ok
  
  list(
    ref  = ref2[keep, , drop = FALSE],
    test = tst2[keep, , drop = FALSE]
  )
}

evalue_rna <- function(path_real, path_impute, meta_path, group,
                       name_real, name_impute, name_method, feat_reduce = NULL){
  
  stopifnot(length(name_impute) == length(name_method))
  
  # 读入真值（可按需降特征）
  if (!is.null(feat_reduce)) {
    feat <- fread(feat_reduce, data.table = FALSE, header = TRUE)[, 1]
    gt_counts <- as.matrix(readRDS(file.path(path_real, name_real)))[feat, ]
  } else {
    gt_counts <- as.matrix(readRDS(file.path(path_real, name_real)))
  }
  
  # sMSE参照
  r_data <- as.matrix(NormalizeData(gt_counts))
  
  # meta / 细胞类型
  meta <- readRDS(meta_path)
  cty  <- meta[, group]
  
  # 仅依赖真值的一次性准备
  sce <- SingleCellExperiment(list(counts = gt_counts))
  sce <- logNormCounts(sce)
  gt_log <- logcounts(sce)
  
  # pDES: 每类细胞 top5 基因
  limma.res <- doLimma(as.matrix(gt_log), cty)
  gene_list <- lapply(names(limma.res), function(nm) rownames(limma.res[[nm]])[1:5])
  names(gene_list) <- names(limma.res)
  hvgs_des <- Reduce(union, gene_list)
  
  # pFCS: HVG100
  var.out  <- modelGeneVar(sce)
  hvgs_fcs <- getTopHVGs(var.out, n = 100)
  
  # 对齐工具
  align_by_genes <- function(ref_mat, test_mat, genes){
    g <- intersect(genes, intersect(rownames(ref_mat), rownames(test_mat)))
    list(
      ref  = ref_mat[g, , drop = FALSE],
      test = test_mat[g, , drop = FALSE],
      genes = g
    )
  }
  
  # 结果表
  out <- data.frame(sMSE = rep(NA_real_, length(name_method)),
                    pFCS = rep(NA_real_, length(name_method)),
                    pDES = rep(NA_real_, length(name_method)))
  rownames(out) <- name_method
  
  # 主循环
  for (m in seq_along(name_method)) {
    cat("Evaluating:", name_method[m], "\n")
    
    # 读预测
    fp <- file.path(path_impute, name_impute[m])
    pd <- as.matrix(fread(fp))
    
    if (name_method[m] == "Palette") {
      p_data <- pd
    } else if (name_method[m] == "scVAEIT") {
      p_data <- t(pd[-1, , drop = FALSE])
    } else if (name_method[m] == "Multigrate") {
      p_data <- as.matrix(NormalizeData(t(pd[-1, , drop = FALSE])))
    } else {
      p_data <- as.matrix(NormalizeData(t(pd[, -1, drop = FALSE])))
    }
    
    # 对齐维度名（若预测没带 dimnames，就强制对齐真值）
    rownames(p_data) <- rownames(gt_counts)
    colnames(p_data) <- colnames(gt_counts)
    
    # sMSE
    out[m, "sMSE"] <- mean(compute_smse(p_data, r_data), na.rm = TRUE)
    
    ## ===== pDES：清洗后再算相关 =====
    al1  <- align_by_genes(gt_log, p_data, hvgs_des)
    flt1 <- drop_zero_or_const_rows(al1$ref, al1$test)
    
    if (nrow(flt1$ref) >= 2) {
      cor_gt1 <- cor(t(flt1$ref),  use = "pairwise.complete.obs")
      cor_p1  <- cor(t(flt1$test), use = "pairwise.complete.obs")
      
      hc_order <- {
        dd <- stats::as.dist((1 - cor_gt1) / 2)
        stats::hclust(dd, method = "complete")$order
      }
      cor_gt1 <- cor_gt1[hc_order, hc_order]
      cor_p1  <- cor_p1 [hc_order, hc_order]
      
      cor_summary_des <- vapply(seq_len(nrow(cor_p1)), function(j){
        suppressWarnings(cor(cor_gt1[j, ], cor_p1[j, ], use = "pairwise.complete.obs"))
      }, numeric(1))
      out[m, "pDES"] <- mean(cor_summary_des, na.rm = TRUE)
    } else {
      warning(sprintf("Method %s (pDES): too few features after filtering; set NA.", name_method[m]))
      out[m, "pDES"] <- NA_real_
    }
    
    ## ===== pFCS：清洗后再算相关 =====
    al2  <- align_by_genes(gt_log, p_data, hvgs_fcs)
    flt2 <- drop_zero_or_const_rows(al2$ref, al2$test)
    
    if (nrow(flt2$ref) >= 2) {
      cor_gt2 <- cor(t(flt2$ref),  use = "pairwise.complete.obs")
      cor_p2  <- cor(t(flt2$test), use = "pairwise.complete.obs")
      
      hc_order2 <- {
        dd <- stats::as.dist((1 - cor_gt2) / 2)
        stats::hclust(dd, method = "complete")$order
      }
      cor_gt2 <- cor_gt2[hc_order2, hc_order2]
      cor_p2  <- cor_p2 [hc_order2, hc_order2]
      
      cor_summary_fcs <- vapply(seq_len(nrow(cor_p2)), function(j){
        suppressWarnings(cor(cor_gt2[j, ], cor_p2[j, ], use = "pairwise.complete.obs"))
      }, numeric(1))
      out[m, "pFCS"] <- mean(cor_summary_fcs, na.rm = TRUE)
    } else {
      warning(sprintf("Method %s (pFCS): too few features after filtering; set NA.", name_method[m]))
      out[m, "pFCS"] <- NA_real_
    }
  }
  
  out
}

evalue_adt <- function(path_real, path_impute, meta_path, group,
                       name_real, name_impute, name_method){
  
  stopifnot(length(name_impute) == length(name_method))
  
  # ---- 1) 真实 ADT & 两种参照（用于 sMSE 也将用于结构相似） ----
  gt_counts  <- as.matrix(readRDS(file.path(path_real, name_real)))
  r_data_log <- as.matrix(Seurat::NormalizeData(gt_counts))  # LogNormalize (log1p)
  r_data_clr <- as.matrix(Seurat::NormalizeData(gt_counts, normalization.method = "CLR", margin = 2))
  
  # ---- 2) meta/细胞类型 & 仅依赖真值的一次性准备（特征选择用）----
  meta <- readRDS(meta_path)
  cty  <- meta[, group]
  
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = gt_counts))
  sce <- scuttle::logNormCounts(sce)
  gt_log <- SummarizedExperiment::assay(sce, "logcounts")  # 仅用于特征选择
  
  # pDES: 每类细胞 top5 基因（保持不变）
  limma.res <- doLimma(as.matrix(gt_log), cty)
  gene_list <- lapply(names(limma.res), function(nm) rownames(limma.res[[nm]])[1:5])
  names(gene_list) <- names(limma.res)
  hvgs_des <- Reduce(union, gene_list)
  
  # pFCS: 前100 HVG（保持不变）
  var.out  <- scran::modelGeneVar(sce)
  hvgs_fcs <- scran::getTopHVGs(var.out, n = 100)
  
  # ---- 3) 工具 ----
  align_by_genes <- function(ref_mat, test_mat, genes){
    g <- intersect(genes, intersect(rownames(ref_mat), rownames(test_mat)))
    list(
      ref  = ref_mat[g, , drop = FALSE],
      test = test_mat[g, , drop = FALSE],
      genes = g
    )
  }
  
  read_pred_adt <- function(fp, method){
    pd <- data.table::fread(fp); pd <- as.matrix(pd)
    if (method %in% c("Palette")) {
      p_data <- pd
    } else if (method %in% c("scVAEIT", "Multigrate")) {
      p_data <- t(pd[-1, , drop = FALSE])
    } else {
      p_data <- t(pd[, -1, drop = FALSE])
    }
    storage.mode(p_data) <- "double"
    p_data
  }
  
  normalize_to_scale <- function(mat, scale = c("CLR","log1p","none")){
    # scale <- match.arg(scale)
    if (scale == "CLR") {
      as.matrix(Seurat::NormalizeData(mat, normalization.method = "CLR", margin = 2))
    } else if (scale == "log1p") {
      as.matrix(Seurat::NormalizeData(mat))
    } else {
      mat
    }
  }
  
  # 方法 → 目标尺度（既用于 sMSE 参照，也用于 pDES/pFCS 的真实参照）
  method_norm_map <- c(
    Palette    = 'none',
    MIDAS      = "CLR",
    Multigrate = 'none',
    scVAEIT    = "log1p"
  )
  get_target_scale <- function(method){
    if (!is.null(method_norm_map[method])) method_norm_map[method] else "CLR"
  }
  
  # ---- 4) 结果表 ----
  out <- data.frame(sMSE = NA_real_, pFCS = NA_real_, pDES = NA_real_)
  out <- out[rep(1, length(name_method)), , drop = FALSE]
  rownames(out) <- name_method
  
  # ---- 5) 主循环 ----
  for (m in seq_along(name_method)) {
    method <- name_method[m]
    cat("Evaluating:", method, "\n")
    
    # 读入 & 对齐维度名
    fp <- file.path(path_impute, name_impute[m]) 
    p_raw <- read_pred_adt(fp, method)
    
    rownames(p_raw) <- rownames(gt_counts)
    colnames(p_raw) <- colnames(gt_counts)
    
    genes_inter <- intersect(rownames(gt_counts), rownames(p_raw))
    cells_inter <- intersect(colnames(gt_counts), colnames(p_raw))
    gt_counts2 <- gt_counts[genes_inter, cells_inter, drop = FALSE]
    p_raw      <- p_raw     [genes_inter, cells_inter, drop = FALSE]
    print(dim(p_raw))
    # 归一化到方法的目标尺度
    target_scale <- get_target_scale(method)
    p_data <- normalize_to_scale(p_raw, target_scale)
    
    # ---- sMSE：用对应尺度的真实参照 ----
    r_ref <- if (target_scale == "log1p") r_data_log[genes_inter, cells_inter, drop = FALSE]
    else                         r_data_clr[genes_inter, cells_inter, drop = FALSE]
    out[m, "sMSE"] <- mean(compute_smse(p_data, r_ref), na.rm = TRUE)
    
    # ---- ★ pDES/pFCS：相似矩阵也用对应尺度的真实参照 ★ ----
    gt_struct_ref <- if (target_scale == "log1p") r_data_log[genes_inter, cells_inter, drop = FALSE]
    else                         r_data_clr[genes_inter, cells_inter, drop = FALSE]
    
    # pDES（特征保持用 gt_log 选出的 hvgs_des；相似矩阵用 gt_struct_ref）
    al1 <- align_by_genes(gt_struct_ref, p_data, hvgs_des)
    cor_gt1 <- cor(t(al1$ref),  use = "na.or.complete")
    cor_p1  <- cor(t(al1$test), use = "na.or.complete")
    hc_order <- {
      dd <- stats::as.dist((1 - cor_gt1) / 2)
      stats::hclust(dd, method = "complete")$order
    }
    cor_gt1 <- cor_gt1[hc_order, hc_order]
    cor_p1  <- cor_p1 [hc_order, hc_order]
    cor_summary_des <- vapply(seq_len(nrow(cor_p1)), function(j){
      suppressWarnings(cor(cor_gt1[j, ], cor_p1[j, ]))
    }, numeric(1))
    out[m, "pDES"] <- mean(cor_summary_des, na.rm = TRUE)
    
    # pFCS（特征保持 hvgs_fcs；相似矩阵用同一个 gt_struct_ref）
    al2 <- align_by_genes(gt_struct_ref, p_data, hvgs_fcs)
    cor_gt2 <- cor(t(al2$ref),  use = "na.or.complete")
    cor_p2  <- cor(t(al2$test), use = "na.or.complete")
    hc_order2 <- {
      dd <- stats::as.dist((1 - cor_gt2) / 2)
      stats::hclust(dd, method = "complete")$order
    }
    cor_gt2 <- cor_gt2[hc_order2, hc_order2]
    cor_p2  <- cor_p2 [hc_order2, hc_order2]
    cor_summary_fcs <- vapply(seq_len(nrow(cor_p2)), function(j){
      suppressWarnings(cor(cor_gt2[j, ], cor_p2[j, ]))
    }, numeric(1))
    out[m, "pFCS"] <- mean(cor_summary_fcs, na.rm = TRUE)
  }
  
  out
}

evalue_atac <- function(path_real, path_impute, meta_path, group,
                        name_real, name_impute, name_method, feat_reduce = NULL){
  
  stopifnot(length(name_impute) == length(name_method))
  if(!is.null(feat_reduce)){
    feat <- fread(feat_reduce,data.table = F, header = T)[,1]
    gt_counts  <- as(as.matrix(readRDS(file.path(path_real, name_real))),'dgCMatrix')[feat,]
  }else{
    gt_counts  <- as(as.matrix(readRDS(file.path(path_real, name_real))),'dgCMatrix')
  }
  # ---- 1) 真实 ATAC & 两种参照（用于 sMSE 也将用于结构相似） ----
  # gt_counts  <- as(as.matrix(readRDS(file.path(path_real, name_real))),'dgCMatrix')
  # r_data_log <- Seurat::NormalizeData(gt_counts)
  
  # Palette
  r_data_idf <- IDFLog(gt_counts)
  
  # MIDAS and Multigrate
  r_data_log <- Seurat::NormalizeData(gt_counts)
  
  # scVAEIT
  r_data_bi <- gt_counts
  r_data_bi@x <- rep(1,length(gt_counts@x))
  
  # ---- 2) meta/细胞类型 & 仅依赖真值的一次性准备（特征选择用）----
  meta <- readRDS(meta_path)
  cty  <- meta[, group]
  
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = gt_counts))
  sce <- scuttle::logNormCounts(sce)
  gt_log <- SummarizedExperiment::assay(sce, "logcounts")  # 仅用于特征选择
  
  # pDES: 每类细胞 top5 基因（保持不变）
  limma.res <- doLimma(as.matrix(gt_log), cty)
  gene_list <- lapply(names(limma.res), function(nm) rownames(limma.res[[nm]])[1:50])
  names(gene_list) <- names(limma.res)
  hvgs_des <- Reduce(union, gene_list)
  
  # pFCS: 前100 HVG（保持不变）
  var.out  <- scran::modelGeneVar(sce)
  hvgs_fcs <- scran::getTopHVGs(var.out, n = 1000)
  
  # ---- 3) 工具 ----
  align_by_genes <- function(ref_mat, test_mat, genes){
    g <- intersect(genes, intersect(rownames(ref_mat), rownames(test_mat)))
    list(
      ref  = ref_mat[g, , drop = FALSE],
      test = test_mat[g, , drop = FALSE],
      genes = g
    )
  }
  
  read_pred_adt <- function(fp, method){
    if(method %in% c("Palette",'Multigrate')){
      pd <- Matrix::readMM(fp)
    }else{
      pd <- data.table::fread(fp)
      pd <- as(as.matrix(pd),'dgCMatrix')
    }
    
    if (method %in% c("Palette")) {
      p_data <- pd
    } else if (method %in% c("scVAEIT")) {
      p_data <- t(pd[-1, , drop = FALSE])
    } else if(method %in% c("Multigrate")){
      p_data <- t(pd)
    }else {
      p_data <- t(pd[, -1, drop = FALSE])
    }
    # storage.mode(p_data) <- "double"
    p_data
  }
  
  normalize_to_scale <- function(mat, scale = c("log1p","none")){
    if (scale == "log1p") {
      as.matrix(Seurat::NormalizeData(mat))
    } else {
      mat
    }
  }
  
  # 方法 → 目标尺度（既用于 sMSE 参照，也用于 pDES/pFCS 的真实参照）
  method_norm_map <- c(
    Palette    = 'none',
    MIDAS      = "log1p",
    Multigrate = 'none',
    scVAEIT    = "none"
  )
  get_target_scale <- function(method){
    if (!is.null(method_norm_map[method])) method_norm_map[method] else "none"
  }
  
  # ---- 4) 结果表 ----
  out <- data.frame(sMSE = NA_real_, pFCS = NA_real_, pDES = NA_real_)
  out <- out[rep(1, length(name_method)), , drop = FALSE]
  rownames(out) <- name_method
  
  # ---- 5) 主循环 ----
  for (m in seq_along(name_method)) {
    method <- name_method[m]
    cat("Evaluating:", method, "\n")
    
    # 读入 & 对齐维度名
    fp <- file.path(path_impute, name_impute[m]) 
    p_raw <- read_pred_adt(fp, method)
    
    rownames(p_raw) <- rownames(gt_counts)
    colnames(p_raw) <- colnames(gt_counts)
    
    genes_inter <- intersect(rownames(gt_counts), rownames(p_raw))
    cells_inter <- intersect(colnames(gt_counts), colnames(p_raw))
    gt_counts2 <- gt_counts[genes_inter, cells_inter, drop = FALSE]
    p_raw      <- p_raw[genes_inter, cells_inter, drop = FALSE]
    print(dim(p_raw))
    # 归一化到方法的目标尺度
    target_scale <- get_target_scale(method)
    p_data <- normalize_to_scale(p_raw, target_scale)
    
    # ---- sMSE：用对应尺度的真实参照 ----
    if(method == 'Palette'){
      r_ref <- r_data_idf[genes_inter, cells_inter, drop = FALSE]
    }
    if(method %in% c("MIDAS",'Multigrate')){
      r_ref <- r_data_log[genes_inter, cells_inter, drop = FALSE]
    }
    if(method == 'scVAEIT'){
      r_ref <- r_data_bi[genes_inter, cells_inter, drop = FALSE]
    }
    out[m, "sMSE"] <- mean(compute_smse(as.matrix(p_data), 
                                        as.matrix(r_ref)), na.rm = TRUE)
    
    # ---- ★ pDES/pFCS：相似矩阵也用对应尺度的真实参照 ★ ----
    if(method == 'Palette'){
      gt_struct_ref <- r_data_idf[genes_inter, cells_inter, drop = FALSE]
    }
    if(method %in% c("MIDAS",'Multigrate')){
      gt_struct_ref <- r_data_log[genes_inter, cells_inter, drop = FALSE]
    }
    if(method == 'scVAEIT'){
      gt_struct_ref <- r_data_bi[genes_inter, cells_inter, drop = FALSE]
    }
    
    # pDES（特征保持用 gt_log 选出的 hvgs_des；相似矩阵用 gt_struct_ref）
    al1 <- align_by_genes(gt_struct_ref, p_data, hvgs_des)
    flt1 <- drop_zero_or_const_rows(al1$ref, al1$test)
    
    if (nrow(flt1$ref) >= 2) {
      cor_gt1 <- cor(t(flt1$ref),  use = "pairwise.complete.obs")
      cor_p1  <- cor(t(flt1$test), use = "pairwise.complete.obs")
      
      hc_order <- {
        dd <- stats::as.dist((1 - cor_gt1) / 2)
        stats::hclust(dd, method = "complete")$order
      }
      cor_gt1 <- cor_gt1[hc_order, hc_order]
      cor_p1  <- cor_p1 [hc_order, hc_order]
      
      cor_summary_des <- vapply(seq_len(nrow(cor_p1)), function(j){
        suppressWarnings(cor(cor_gt1[j, ], cor_p1[j, ], use = "pairwise.complete.obs"))
      }, numeric(1))
      out[m, "pDES"] <- mean(cor_summary_des, na.rm = TRUE)
    } else {
      warning(sprintf("Method %s (pDES): too few features after filtering; set NA.", method))
      out[m, "pDES"] <- NA_real_
    }
    
    ############
    # cor_gt1 <- cor(as.matrix(t(al1$ref)),  use = "na.or.complete")
    # cor_p1  <- cor(as.matrix(t(al1$test)), use = "na.or.complete")
    # hc_order <- {
    #   dd <- stats::as.dist((1 - cor_gt1) / 2)
    #   stats::hclust(dd, method = "complete")$order
    # }
    # cor_gt1 <- cor_gt1[hc_order, hc_order]
    # cor_p1  <- cor_p1 [hc_order, hc_order]
    # cor_summary_des <- vapply(seq_len(nrow(cor_p1)), function(j){
    #   suppressWarnings(cor(cor_gt1[j, ], cor_p1[j, ]))
    # }, numeric(1))
    # out[m, "pDES"] <- mean(cor_summary_des, na.rm = TRUE)
    
    
    # pFCS（特征保持 hvgs_fcs；相似矩阵用同一个 gt_struct_ref）
    al2 <- align_by_genes(gt_struct_ref, p_data, hvgs_fcs)
    flt2 <- drop_zero_or_const_rows(al2$ref, al2$test)
    
    if (nrow(flt2$ref) >= 2) {
      cor_gt2 <- cor(t(flt2$ref),  use = "pairwise.complete.obs")
      cor_p2  <- cor(t(flt2$test), use = "pairwise.complete.obs")
      
      hc_order2 <- {
        dd <- stats::as.dist((1 - cor_gt2) / 2)
        stats::hclust(dd, method = "complete")$order
      }
      cor_gt2 <- cor_gt2[hc_order2, hc_order2]
      cor_p2  <- cor_p2 [hc_order2, hc_order2]
      
      cor_summary_fcs <- vapply(seq_len(nrow(cor_p2)), function(j){
        suppressWarnings(cor(cor_gt2[j, ], cor_p2[j, ], use = "pairwise.complete.obs"))
      }, numeric(1))
      out[m, "pFCS"] <- mean(cor_summary_fcs, na.rm = TRUE)
    } else {
      warning(sprintf("Method %s (pFCS): too few features after filtering; set NA.", method))
      out[m, "pFCS"] <- NA_real_
    }
    ###########
  #   cor_gt2 <- cor(as.matrix(t(al2$ref)),  use = "na.or.complete")
  #   cor_p2  <- cor(as.matrix(t(al2$test)), use = "na.or.complete")
  #   hc_order2 <- {
  #     dd <- stats::as.dist((1 - cor_gt2) / 2)
  #     stats::hclust(dd, method = "complete")$order
  #   }
  #   cor_gt2 <- cor_gt2[hc_order2, hc_order2]
  #   cor_p2  <- cor_p2 [hc_order2, hc_order2]
  #   cor_summary_fcs <- vapply(seq_len(nrow(cor_p2)), function(j){
  #     suppressWarnings(cor(cor_gt2[j, ], cor_p2[j, ]))
  #   }, numeric(1))
  #   out[m, "pFCS"] <- mean(cor_summary_fcs, na.rm = TRUE)
  }
  
  out
}
IDFLog <- function(data,
                   scale.factor = 1e4,
                   verbose = TRUE){
  rsums <- rowSums(x = data)
  col_sum <- ncol(data)
  idf <- col_sum / rsums
  norm.data <- Matrix::Diagonal(n = length(x = idf), x = idf)%*%data
  norm.data <- Seurat::NormalizeData(norm.data,
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
  return(norm.data)
}
