library(StabMap)
library(scran)
library(Seurat)
library(patchwork)

setwd('./10x_Visium_human_tonsil')
out_dir <- './10xVisium_tonsi'

paths <- c('slice1','slice2','slice3')
adt.list <- list()
rna.list <- list()
meta.list <- list()
for(i in 1:3){
  if(i == 1){
    rna.list[[i]] <- readRDS(paste0('./',paths[i],'/rna.rds'))   
    adt.list[[i]] <- readRDS(paste0('./',paths[i],'/adt.rds'))
  }
  if(i == 2){
    rna.list[[i]] <- readRDS(paste0('./',paths[i],'/rna.rds'))
  }
  if(i == 3 ){
    adt.list[[i-1]] <- readRDS(paste0('./',paths[i],'/adt.rds'))
  }
  meta.list[[i]] <- readRDS(paste0('./',paths[i],'/meta.rds'))
}

for(i in 1:2){
  rna.list[[i]] <- NormalizeData(rna.list[[i]])
  adt.list[[i]] <- NormalizeData(adt.list[[i]])
}

# ADT
adt_a1 <- adt.list[[1]]
adt_a3 <- adt.list[[2]]

rna_a1 <- rna.list[[1]]
rna_a2 <- rna.list[[2]]

rownames(adt_a1) <- paste0(rownames(adt_a1),'-adt')
rownames(adt_a3) <- paste0(rownames(adt_a3),'-adt')

rownames(rna_a1) <- paste0(rownames(rna_a1),'-rna')
rownames(rna_a2) <- paste0(rownames(rna_a2),'-rna')

assay_list = list(B1=rbind(rna_a1,adt_a1),B2=rna_a2,B3=adt_a3)
mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("B1"),
               maxFeatures = 10000,
               ncomponentsReference = 20,
               ncomponentsSubset = 20)

meta <- Reduce(rbind,meta.list)
stab <- stab[rownames(meta),]
write.csv(as.matrix(stab),
          paste0(out_dir,'/stab.csv'),
          row.names = FALSE)
