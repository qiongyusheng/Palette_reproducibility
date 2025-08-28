library(StabMap)
library(scran)
library(Seurat)
library(patchwork)

setwd(".../")
rna.list <- list()
atac.list <- list()
adt.list <- list()
meta.list <- list()
for(i in 2:5){
  rna.list[[i-1]] <- readRDS(paste0("./B",i,"/rna.rds"))
  meta.list[[i-1]] <- readRDS(paste0("./B",i,"/meta.rds"))
  atac.list[[i-1]] <- readRDS(paste0("./B",i,"/atac.rds"))
  adt.list[[i-1]] <- readRDS(paste0("./B",i,"/adt.rds"))
}

atac.list <- atac.list[1:2]
rna.list <- rna.list[c(1,3)]
adt.list <- adt.list[c(1,4)]

for(i in 1:2){
  atac.list[[i]] <- NormalizeData(atac.list[[i]])
  rna.list[[i]] <- NormalizeData(rna.list[[i]])
  adt.list[[i]] <- NormalizeData(adt.list[[i]])
}

# ATAC
d1 <- atac.list[[1]]
d2 <- atac.list[[2]]
# RNA
d3 <- rna.list[[1]]
d4 <- rna.list[[2]]
# ADT
d5 <- adt.list[[1]]
d6 <- adt.list[[2]]

rm(atac.list,rna.list,adt.list )

rownames(d1) <- paste0(rownames(d1),'-atac')
rownames(d2) <- paste0(rownames(d2),'-atac')
rownames(d3) <- paste0(rownames(d3),'-rna')
rownames(d4) <- paste0(rownames(d4),'-rna')
rownames(d5) <- paste0(rownames(d5),'-adt')
rownames(d6) <- paste0(rownames(d6),'-adt')

b1 <- Reduce(rbind,list(d1,d3,d5))
rm(d1,d3,d5)

assay_list = list(Batch1=b1,Batch2=d2,Batch3=d4,Batch4=d6)
#mf <- nrow(d1)
rm(d2,d3,d4,d5,d6,b1)
gc()
mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("Batch1"),
               maxFeatures = nrow(assay_list[[1]]),
               ncomponentsReference = 20,
               ncomponentsSubset = 20)
meta <- readRDS('./meta.rds')
write.csv(as.matrix(stab)[,rownames(meta)],
          "./stab.csv",
          row.names = FALSE)