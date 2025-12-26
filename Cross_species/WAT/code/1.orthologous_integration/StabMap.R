library(StabMap)
library(scran)
library(Seurat)
library(patchwork)

setwd('./human_mouse_WAT')
out_dir <- './'
meta <- readRDS('./meta.rds')
human_homo <- as(as.matrix(readRDS('./human/rna_homo.rds')),'dgCMatrix')
# human_non <- as(as.matrix(readRDS('./human/rna_non.rds')),'dgCMatrix')
# human_gene <- readRDS('./human/gene_non.rds')

mouse_homo <- as(as.matrix(readRDS('./mouse/rna_homo.rds')),'dgCMatrix')
# mouse_non <- as(as.matrix(readRDS('./mouse/rna_non.rds')),'dgCMatrix')
# mouse_gene <- readRDS('./mouse/gene_non.rds')

human <- human_homo
rm(human_homo)
mouse <- mouse_homo
rm(mouse_homo)
human_meta <- readRDS('./human/meta.rds')
mouse_meta <- readRDS('./mouse/meta.rds')
human <- NormalizeData(human)
mouse <- NormalizeData(mouse)

human_list <- list()
for(i in unique(human_meta$Batch)){
    sub_meta <- dplyr::filter(human_meta,Batch == i)
    human_list[[i]] <- human[,rownames(sub_meta)]
}

mouse_list <- list()
for(i in unique(mouse_meta$Batch)){
    sub_meta <- dplyr::filter(mouse_meta,Batch == i)
    mouse_list[[i]] <- mouse[,rownames(sub_meta)]
}

dim(human)
dim(mouse)

rm(human,mouse)

table(human_meta$Batch)
table(mouse_meta$Batch)

assay_list = c(human_list,mouse_list)
mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("Hs254"),
               maxFeatures = 10000,
              ncomponentsReference = 20,
              ncomponentsSubset = 20)

stab <- stab[rownames(meta),]

write.csv(as.matrix(stab),
          paste0(out_dir,'/stab_homo.csv'),
          row.names = FALSE)