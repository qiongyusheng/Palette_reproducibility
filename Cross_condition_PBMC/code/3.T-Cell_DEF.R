library(future)
library(Seurat)
library(data.table)
library(dplyr)
library(purrr)
library(Signac)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(Matrix)
library(matrixStats)
library(edgeR)

setwd("./")
path <- c('./PBMC_ASAP/control_asap','./PBMC_ASAP/stim_asap',
          './PBMC_CITE/control_cite','./PBMC_CITE/stim_cite')
adt.list <- list()
meta.list <- list()
for(i in 1:4){
    meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
    adt.list[[i]] <- readRDS(paste0(path[i],'/adt.rds'))
}

mm <- data.frame(batch = c(meta.list[[1]]$batch,meta.list[[2]]$batch,
                           meta.list[[3]]$batch,meta.list[[4]]$batch),
                 modality = c(meta.list[[1]]$modality,meta.list[[2]]$modality,
                           meta.list[[3]]$modality,meta.list[[4]]$modality),
                 condition = c(meta.list[[1]]$condition,meta.list[[2]]$condition,
                           meta.list[[3]]$condition,meta.list[[4]]$condition),
                  celltype = c(as.character(meta.list[[1]][,1]),
                               as.character(meta.list[[2]][,1]),
                               as.character(meta.list[[3]][,1]),
                              as.character(meta.list[[4]][,1]))
                  )
obj <- CreateSeuratObject(Reduce(cbind,adt.list))
obj[['batch']] <- mm[['batch']]
obj[['modality']] <- mm[['modality']]
obj[['condition']] <- mm[['condition']]
obj[['celltype']] <- mm[['celltype']]
rownames(obj@meta.data) <- colnames(Reduce(cbind,adt.list))

setwd('./Cross_condition_PBMC/output/Embedding')
int <- as.matrix(fread('./Palette.csv',data.table = F))
rownames(int) <- colnames(Reduce(cbind,adt.list))
colnames(int) <- paste0("dim_",seq(ncol(int)))
obj[["Palette"]] <- CreateDimReducObject(embeddings = int, 
                                         key = "Palette_")

obj <- RunUMAP(obj, reduction = 'Palette', dims = 1:20, reduction.name = 'Palette.umap')
obj <- FindNeighbors(obj, 
                     reduction = 'Palette',
                     dims = 1:20) %>% FindClusters(., 
                                                   resolution = 0.8,cluster.name = "Palette")

mm$cluster <- obj@meta.data$RNA_snn_res.0.8
rownames(mm) <- colnames(Reduce(cbind,adt.list))

mm_adt_c <- rownames(dplyr::filter(mm, cluster %in% c(13,0,14,1,4,12,16,2)))
mm_adt_s <- rownames(dplyr::filter(mm, cluster %in% c(3,8,10)))
adt_c <- t(Reduce(cbind,adt.list)[,mm_adt_c])
adt_s <- t(Reduce(cbind,adt.list)[,mm_adt_s])

lapply(1:dim(adt_s)[2], function(i){
  wt <- wilcox.test(adt_c[,i], adt_s[,i])
  fc <- log2(mean(adt_s[,i])/mean(adt_c[,i]))
  data.frame(
    marker = colnames(adt_s)[i],
    foldchange = fc,
    pvalue = wt$p.value
  )
}) %>% rbindlist() %>% data.frame() -> adt_df

adt_df$padj <- p.adjust(adt_df$pvalue)

PROTEIN_hits <- adt_df$padj < 0.01 & (abs(adt_df$foldchange) > 0.5)
adt_dge <- adt_df[,1][PROTEIN_hits]


setwd("./")
rna1 <- Read10X("./rna-seq/ctrl")
dge_rna <- readRDS('./DGE_RNA_Tcell_Activation.rds')

colnames(rna1) <- paste0(colnames(rna1),'cite_control')
rna1 <- rna1[ ,colnames(adt.list[[3]])]

rna2 <- Read10X("./rna-seq/stim")
colnames(rna2) <- paste0(colnames(rna2),'cite_stim')
rna2 <- rna2[ ,colnames(adt.list[[4]])]

mm_rna_c <- rownames(dplyr::filter(dplyr::filter(mm,cluster %in% c(13,0,14,1,4,12,16,2)),modality == 'cite'))
mm_rna_s <- rownames(dplyr::filter(dplyr::filter(mm,cluster %in% c(3,8,10)),modality == 'cite'))

run_edgeRQLFdetect_CL <- function(count, condt) {
  
  dge <- DGEList(count, group = condt)
  dge <- calcNormFactors(dge)
  
  # adjust for total expression
  cdr <- scale(colMeans(count > 0))
  cdr2 <- scale(colSums(count))
  
  design <- model.matrix(~ cdr + condt) 
  dge <- estimateDisp(dge, design = design)
  fit <- glmQLFit(dge, design = design)
  qlf <- glmQLFTest(fit)
  tt <- topTags(qlf, n = Inf, sort.by = "none")
  df <- signif(tt@.Data[[1]], 3)
  df$gene <- rownames(df)
  
  # small function to pull CPM
  ps <- function(which_ones){
    rs <- rowSums(count[,which_ones])
    cpms <- round(rs/sum(rs) *1000000,1)[as.character(df$gene)]
    return(cpms)
  } 
  
  # Pull CPM values
  df$control_cpm <- ps(condt == "c")
  df$stim_cpm <- ps(condt == "s")
  
  # Round
  df$logFC <- round(df$logFC,2)
  df$logCPM <- round(df$logCPM,2)
  df
}

condt <- c(
  rep("c", dim(t(cbind(rna1,rna2)[,mm_rna_c]))[1]),
  rep("s", dim(t(cbind(rna1,rna2)[,mm_rna_s]))[1])
)

counts <- t(data.matrix(rbind(t(cbind(rna1,rna2)[,mm_rna_c]), t(cbind(rna1,rna2)[,mm_rna_s]))))
rs <- rowSums(counts)
cpms <- round(rs/sum(rs) *1000000,1)
counts2 <- counts[cpms > 2,]
RNA_DGE <- run_edgeRQLFdetect_CL(counts2, condt)
RNA_DGE_rank <- RNA_DGE %>% arrange(desc(F))

setwd("./")
atac1 <- as(as.matrix(readRDS('./PBMC_ASAP/control_asap/atac.rds')),'dgCMatrix')
atac2 <- as(as.matrix(readRDS('./PBMC_ASAP/stim_asap/atac.rds')),'dgCMatrix')

mm_atac_c <- rownames(dplyr::filter(dplyr::filter(mm,cluster %in% c(13,0,14,1,4,12,16,2)),modality == 'asap'))
mm_atac_s <- rownames(dplyr::filter(dplyr::filter(mm,cluster %in% c(3,8,10)),modality == 'asap'))

atac_control <- cbind(atac1,atac2)[,mm_atac_c ]
atac_stim <- cbind(atac1,atac2)[,mm_atac_s]

cpm <- function(counts){
  rs <- rowSums(counts)
  cpms <- round(rs/sum(rs) *1000000,1)
  cpms
}
cpms <- cpm(cbind(atac_control,atac_stim))
atac_control <- atac_control[cpms > 2,]
atac_stim <- atac_stim[cpms > 2,]

diff_atac_test <- function(mat1, mat2){
  
  mat <- cbind(mat1, mat2)
  dim(mat)
  
  vec <- c(rep(TRUE, dim(mat1)[2]), rep(FALSE, dim(mat2)[2]))
  
  # Permuted values
  permuted <- sapply(1:100, function(i){
    set.seed(i)
    vec_p <- sample(vec)
    p1 <- rowMeans(mat[,vec_p] != 0)
    p2 <- rowMeans(mat[,!vec_p] != 0)
    p1-p2
  })
  
  # Establish baseline characteristics for differential analyses
  diffDF <- data.frame(
    peak_idx = seq(1:dim(mat)[1]),
    p1 = rowMeans(mat[,vec] != 0),
    p2 = rowMeans(mat[,!vec] != 0),
    p = rowMeans(mat != 0),
    mean_permuted = rowMeans(permuted),
    sd_permuted = rowSds(permuted)
  )
  
  # Formalize test statistics
  diffDF %>% 
    mutate(diffP = p1 - p2) %>%
    mutate(Z = ifelse(sd_permuted > 0, (diffP-mean_permuted)/sd_permuted,0)) %>%
    mutate(pvalue = 2*pnorm(-abs(Z)))%>%
    mutate(log10p = -1*log10(pvalue)) %>%
    mutate(FDR = p.adjust(pvalue)) %>% 
    mutate(log10FDR = -1*log10(FDR)) -> diffDF
  
  diffDF
}

diff_atac <- diff_atac_test(atac_control,atac_stim)
diff_atac$logFC <- log2((diff_atac$p1 + 0.0001)/(diff_atac$p2+ 0.0001))

saveRDS(diff_atac,"./Cross_condition_PBMC/output/DGE/peak_new.rds")
saveRDS(adt_df,'./Cross_condition_PBMC/output/DGE/adt_new.rds')
saveRDS(RNA_DGE_rank,'./Cross_condition_PBMC/output/DGE/rna_new.rds')
