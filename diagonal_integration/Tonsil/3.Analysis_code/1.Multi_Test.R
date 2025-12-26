library(reshape2)
library(RSpectra)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(tidydr)
library(ggpubr)
library(cowplot) 
library(Matrix)
library(rstatix)
library(data.table)

rna <- readMM('./RNA/tonsil_rna_0510.txt')
gene_name <- fread('./RNA/tonsil_rna_0510_names.csv',
                   data.table = F,header = T)
rna_meta <- fread('./RNA/tonsil_rna_0510_meta.csv',
                  data.table = F)
pred_rna <- readRDS('./pred_label_all_methods_RNA.rds')

rna <- as(t(rna),'dgCMatrix')
rownames(rna) <- gene_name[,2]
colnames(rna) <- rna_meta[,1]

setwd('.')
rnfer_data <- readRDS('./infer_rna.rds')
colnames(rnfer_data) <- paste0('cell_',colnames(rnfer_data))
meta_codex <- readRDS('./CODEX/meta.rds')
pred_codex <- readRDS('./pred_label_all_methods_CODEX.rds')

rna <- NormalizeData(rna)
rownames(rna_meta) <- rna_meta[,1]
rna_obj <- CreateSeuratObject(rna,meta.data = rna_meta)
Idents(rna_obj) <- rna_obj[['cluster.info']]
marker_rna <- FindAllMarkers(rna_obj,logfc.threshold = 0.25)
marker_rna <- dplyr::filter(marker_rna,p_val_adj < 0.05)

rna_obj[['pred']] <- pred_rna$Palette
Idents(rna_obj) <- rna_obj[['pred']]
marker_rnap <- FindAllMarkers(rna_obj,logfc.threshold = 0.25)
marker_rnap <- dplyr::filter(marker_rna,p_val_adj < 0.05)


options(future.globals.maxSize = 10 * 1024^3)  # 10 GiB
rownames(meta_codex) <- colnames(rnfer_data)
codex_obj <- CreateSeuratObject(rnfer_data,meta.data = meta_codex)
Idents(codex_obj) <- codex_obj[['celltype']]
marker_codex <- FindAllMarkers(codex_obj,logfc.threshold = 0.25)
marker_codex <- dplyr::filter(marker_codex,p_val_adj < 0.05)

codex_obj[['pred']] <- pred_codex$Palette
Idents(codex_obj) <- codex_obj[['pred']]
marker_codexp <- FindAllMarkers(codex_obj,logfc.threshold = 0.25)
marker_codexp <- dplyr::filter(marker_codex,p_val_adj < 0.05)

marker_rna$group <- paste0(marker_rna$cluster,'_scRNA')
marker_rnap$group <- paste0(marker_rnap$cluster,'_scRNAp')
marker_codex$group <- paste0(marker_codex$cluster,'_codex')
marker_codexp$group <- paste0(marker_codexp$cluster,'_codexp')
res <- rbind(marker_rna,marker_rnap,marker_codex,marker_codexp)
library(SuperExactTest)
combined_marker_sets <- split(res$gene, res$group)


res <- supertest(combined_marker_sets, n=nrow(rna), degree=4)
res <- summary(res)$Table
rownames(res) <- NULL
saveRDS(res,'../multi_test.rds')