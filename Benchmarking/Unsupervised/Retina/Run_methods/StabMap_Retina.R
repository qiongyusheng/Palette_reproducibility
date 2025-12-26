######################### Select Peaks ################################
library(scran)
library(Seurat)
library(patchwork)
library(data.table)
library(Matrix)
library(SingleCellMultiModal)
setwd('/home/server/sqy/data/mosaic/input/Human_Retina')
path <- c('./LGS1OD','./LGS1OS','./LGS2OD','./LGS2OS','./LGS3OD','./LGS3OS','./LVG1OD','./LVG1OS')

atac.list <- list()
for(i in 1:8){
  if(i %in% c(1,2,6,7,8)){
    peak <- fread(paste0(path[i],'/peak.csv'),data.table = F)
    cells <- fread(paste0(path[i],'/bcd.csv'),data.table = F)
    atac.list[[i]] <- readMM(paste0(path[i],'/atac.mtx'))
    colnames(atac.list[[i]]) <- cells[,1]
    rownames(atac.list[[i]]) <- peak[,1]
  }else{
    atac.list[[i]] <- NULL 
  }
}

atac.list <- atac.list[c(1,2,6,7,8)]
all_atac <- Reduce(cbind,atac.list)
all_atac <- SingleCellExperiment(assays = list(counts = all_atac))
all_atac <- logNormCounts(all_atac)
decomp <- modelGeneVar(all_atac)

atac_hvgs <- rownames(decomp)[decomp$mean > 0.005 & decomp$p.value <= 0.05]
saveRDS(decomp,'./stabma_hvg.rds')


######################### Random 1 ################################
library(StabMap)
library(scran)
library(Seurat)
library(patchwork)
library(data.table)
library(Matrix)
library(SingleCellMultiModal)

setwd('/home/server/sqy/data/mosaic/input/Human_Retina')
path <- c('./LGS1OD','./LGS1OS','./LGS2OD','./LGS2OS','./LGS3OD','./LGS3OS','./LVG1OD','./LVG1OS')

rna.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:8){
  if(i %in% c(3,4,5,6,8)){
    rna.list[[i]] <- readRDS(paste0(path[i],'/rna.rds')) 
  }else{
    rna.list[[i]] <- NULL 
  }
  
  if(i %in% c(1,2,4,6,7)){
    peak <- fread(paste0(path[i],'/peak.csv'),data.table = F)
    cells <- fread(paste0(path[i],'/bcd.csv'),data.table = F)
    atac.list[[i]] <- readMM(paste0(path[i],'/atac.mtx'))
    colnames(atac.list[[i]]) <- cells[,1]
    rownames(atac.list[[i]]) <- peak[,1]
  }else{
    atac.list[[i]] <- NULL 
  }
  meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
}

decomp <- readRDS('./stabma_hvg.rds')
atac_hvgs <- rownames(decomp)[decomp$mean > 0.005 & decomp$p.value <= 0.05]

for(i in 1:8){
  
  if(i %in% c(3,4,5,6,8)){
    rna.list[[i]] <- NormalizeData(rna.list[[i]]) 
  }
  if(i %in% c(1,2,4,6,7)){
    atac.list[[i]] <- NormalizeData(atac.list[[i]][atac_hvgs,])
  }
}

assay_list = list(Batch1=atac.list[[1]],
                  Batch2=atac.list[[2]],
                  Batch3=rna.list[[3]],
                  Batch4=rbind(atac.list[[4]],rna.list[[4]]),
                  Batch5=rna.list[[5]],
                  Batch6=rbind(atac.list[[6]],rna.list[[6]]),
                  Batch7=atac.list[[7]],
                  Batch8=rna.list[[8]])
mosaicDataUpSet(assay_list, plot = F)
rm(rna.list,atac.list)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("Batch4"),
               maxFeatures = nrow(assay_list[[4]]),
               ncomponentsReference = 20,
               ncomponentsSubset = 20)

mm <- Reduce(rbind,meta.list)
emb <- as.matrix(stab)
emb <- emb[rownames(mm),]

write.csv(as.matrix(emb),
          "./Benchmarking/Unsupervised/Retina/output/stab_s1.csv",
          row.names = FALSE)

######################### Random 2 ################################
library(StabMap)
library(scran)
library(Seurat)
library(patchwork)
library(data.table)
library(Matrix)
library(SingleCellMultiModal)

setwd('/home/server/sqy/data/mosaic/input/Human_Retina')
path <- c('./LGS1OD','./LGS1OS','./LGS2OD','./LGS2OS','./LGS3OD','./LGS3OS','./LVG1OD','./LVG1OS')

rna.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:8){
  if(i %in% c(2,3,4,5,7)){
    rna.list[[i]] <- readRDS(paste0(path[i],'/rna.rds')) 
  }else{
    rna.list[[i]] <- NULL 
  }
  
  if(i %in% c(1,2,3,6,8)){
    peak <- fread(paste0(path[i],'/peak.csv'),data.table = F)
    cells <- fread(paste0(path[i],'/bcd.csv'),data.table = F)
    atac.list[[i]] <- readMM(paste0(path[i],'/atac.mtx'))
    colnames(atac.list[[i]]) <- cells[,1]
    rownames(atac.list[[i]]) <- peak[,1]
  }else{
    atac.list[[i]] <- NULL 
  }
  meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
}

decomp <- readRDS('./stabma_hvg.rds')
atac_hvgs <- rownames(decomp)[decomp$mean > 0.005 & decomp$p.value <= 0.05]

for(i in 1:8){
  
  if(i %in% c(2,3,4,5,7)){
    rna.list[[i]] <- NormalizeData(rna.list[[i]]) 
  }
  if(i %in% c(1,2,3,6,8)){
    atac.list[[i]] <- NormalizeData(atac.list[[i]][atac_hvgs,])
  }
}

assay_list = list(Batch1=atac.list[[1]],
                  Batch2=rbind(atac.list[[2]],rna.list[[2]]),
                  Batch3=rbind(atac.list[[3]],rna.list[[3]]),
                  Batch4=rna.list[[4]],
                  Batch5=rna.list[[5]],
                  Batch6=atac.list[[6]],
                  Batch7=rna.list[[7]],
                  Batch8=atac.list[[8]])
mosaicDataUpSet(assay_list, plot = F)
rm(rna.list,atac.list)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("Batch2"),
               maxFeatures = nrow(assay_list[[3]]),
               ncomponentsReference = 20,
               ncomponentsSubset = 20)

mm <- Reduce(rbind,meta.list)
emb <- as.matrix(stab)
emb <- emb[rownames(mm),]

write.csv(as.matrix(emb),
          "./Benchmarking/Unsupervised/Retina/output/stab_s2.csv",
          row.names = FALSE)

######################### Random 3 ################################
library(StabMap)
library(scran)
library(Seurat)
library(patchwork)
library(data.table)
library(Matrix)
library(SingleCellMultiModal)

setwd('/home/server/sqy/data/mosaic/input/Human_Retina')
path <- c('./LGS1OD','./LGS1OS','./LGS2OD','./LGS2OS','./LGS3OD','./LGS3OS','./LVG1OD','./LVG1OS')

rna.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:8){
  if(i %in% c(1,2,3,4,5)){
    rna.list[[i]] <- readRDS(paste0(path[i],'/rna.rds')) 
  }else{
    rna.list[[i]] <- NULL 
  }
  
  if(i %in% c(1,2,6,7,8)){
    peak <- fread(paste0(path[i],'/peak.csv'),data.table = F)
    cells <- fread(paste0(path[i],'/bcd.csv'),data.table = F)
    atac.list[[i]] <- readMM(paste0(path[i],'/atac.mtx'))
    colnames(atac.list[[i]]) <- cells[,1]
    rownames(atac.list[[i]]) <- peak[,1]
  }else{
    atac.list[[i]] <- NULL 
  }
  meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
}

rna.list <- rna.list[c(c(1,2,3,4,5))]
atac.list <- atac.list[c(1,2,6,7,8)]

decomp <- readRDS('./stabma_hvg.rds')
atac_hvgs <- rownames(decomp)[decomp$mean > 0.005 & decomp$p.value <= 0.05]

for(i in 1:5){
  atac.list[[i]] <- NormalizeData(atac.list[[i]][atac_hvgs,])
  rna.list[[i]] <- NormalizeData(rna.list[[i]])
}

assay_list = list(Batch1=rbind(rna.list[[1]],atac.list[[1]]),
                  Batch2=rbind(rna.list[[2]],atac.list[[2]]),
                  Batch3=rna.list[[3]],
                  Batch4=rna.list[[4]],
                  Batch5=rna.list[[5]],
                  Batch6=atac.list[[3]],
                  Batch7=atac.list[[4]],
                  Batch8=atac.list[[5]])
mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("Batch2"),
               maxFeatures = nrow(assay_list[[1]]),
               ncomponentsReference = 20,
               ncomponentsSubset = 20)

mm <- data.frame(batch = c(rep("B1",8455),rep("B2",9182),rep("B3",6234),rep("B4",5955),
                           rep("B5",5707),rep("B6",5013),rep("B7",5037),rep("B8",4729)),
                 modality = c(rep('multiome',8455+9182),rep('rna',6234+5955+5707),rep('atac',5013+5037+4729)),
                 celltype = c(as.character(meta.list[[1]][,4]),
                              as.character(meta.list[[2]][,4]),
                              as.character(meta.list[[3]][,4]),
                              as.character(meta.list[[4]][,4]),
                              as.character(meta.list[[5]][,4]),
                              as.character(meta.list[[6]][,4]),
                              as.character(meta.list[[7]][,4]),
                              as.character(meta.list[[8]][,4]))
)

write.csv(as.matrix(stab),
          "./Benchmarking/Unsupervised/Retina/output/stab.csv",
          row.names = FALSE)

######################### Random 4 ################################
library(StabMap)
library(scran)
library(Seurat)
library(patchwork)
library(data.table)
library(Matrix)
library(SingleCellMultiModal)

setwd('/home/server/sqy/data/mosaic/input/Human_Retina')
path <- c('./LGS1OD','./LGS1OS','./LGS2OD','./LGS2OS','./LGS3OD','./LGS3OS','./LVG1OD','./LVG1OS')

rna.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:8){
  if(i %in% c(2,4,5,6,7)){
    rna.list[[i]] <- readRDS(paste0(path[i],'/rna.rds')) 
  }else{
    rna.list[[i]] <- NULL 
  }
  
  if(i %in% c(1,3,4,5,8)){
    peak <- fread(paste0(path[i],'/peak.csv'),data.table = F)
    cells <- fread(paste0(path[i],'/bcd.csv'),data.table = F)
    atac.list[[i]] <- readMM(paste0(path[i],'/atac.mtx'))
    colnames(atac.list[[i]]) <- cells[,1]
    rownames(atac.list[[i]]) <- peak[,1]
  }else{
    atac.list[[i]] <- NULL 
  }

  meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
}

decomp <- readRDS('./stabma_hvg.rds')
atac_hvgs <- rownames(decomp)[decomp$mean > 0.005 & decomp$p.value <= 0.05]

for(i in 1:8){
  
  if(i %in% c(2,4,5,6,7)){
    rna.list[[i]] <- NormalizeData(rna.list[[i]]) 
  }
  if(i %in% c(1,3,4,5,8)){
    atac.list[[i]] <- NormalizeData(atac.list[[i]][atac_hvgs,])
  }

}

assay_list = list(Batch1=atac.list[[1]],
                  Batch2=rna.list[[2]],
                  Batch3=atac.list[[3]],
                  Batch4=rbind(rna.list[[4]],atac.list[[4]]),
                  Batch5=rbind(rna.list[[5]],atac.list[[5]]),
                  Batch6=rna.list[[6]],
                  Batch7=rna.list[[7]],
                  Batch8=atac.list[[8]])
mosaicDataUpSet(assay_list, plot = F)
rm(rna.list,atac.list)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("Batch4"),
               maxFeatures = nrow(assay_list[[4]]),
               ncomponentsReference = 20,
               ncomponentsSubset = 20)

mm <- Reduce(rbind,meta.list)
emb <- as.matrix(stab)
emb <- emb[rownames(mm),]

write.csv(as.matrix(emb),
          "./Benchmarking/Unsupervised/Retina/output/stab_s3.csv",
          row.names = FALSE)

######################### Random 5 ################################
library(StabMap)
library(scran)
library(Seurat)
library(patchwork)
library(data.table)
library(Matrix)
library(SingleCellMultiModal)

setwd('/home/server/sqy/data/mosaic/input/Human_Retina')
path <- c('./LGS1OD','./LGS1OS','./LGS2OD','./LGS2OS','./LGS3OD','./LGS3OS','./LVG1OD','./LVG1OS')

rna.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:8){
  if(i %in% c(1,3,5,6,7)){
    rna.list[[i]] <- readRDS(paste0(path[i],'/rna.rds')) 
  }else{
    rna.list[[i]] <- NULL 
  }
  
  if(i %in% c(2,4,6,7,8)){
    peak <- fread(paste0(path[i],'/peak.csv'),data.table = F)
    cells <- fread(paste0(path[i],'/bcd.csv'),data.table = F)
    atac.list[[i]] <- readMM(paste0(path[i],'/atac.mtx'))
    colnames(atac.list[[i]]) <- cells[,1]
    rownames(atac.list[[i]]) <- peak[,1]
  }else{
    atac.list[[i]] <- NULL 
  }

  meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
}

decomp <- readRDS('./stabma_hvg.rds')
atac_hvgs <- rownames(decomp)[decomp$mean > 0.005 & decomp$p.value <= 0.05]


for(i in 1:8){
  
  if(i %in% c(1,3,5,6,7)){
    rna.list[[i]] <- NormalizeData(rna.list[[i]]) 
  }
  if(i %in% c(2,4,6,7,8)){
    atac.list[[i]] <- NormalizeData(atac.list[[i]][atac_hvgs,])
  }
}

assay_list = list(Batch1=rna.list[[1]],
                  Batch2=atac.list[[2]],
                  Batch3=rna.list[[3]],
                  Batch4=atac.list[[4]],
                  Batch5=rna.list[[5]],
                  Batch6=rbind(rna.list[[6]],atac.list[[6]]),
                  Batch7=rbind(rna.list[[7]],atac.list[[7]]),
                  Batch8=atac.list[[8]])
mosaicDataUpSet(assay_list, plot = F)
rm(rna.list,atac.list)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("Batch7"),
               maxFeatures = nrow(assay_list[[7]]),
               ncomponentsReference = 20,
               ncomponentsSubset = 20)

mm <- Reduce(rbind,meta.list)
emb <- as.matrix(stab)
emb <- emb[rownames(mm),]

write.csv(as.matrix(emb),
          "./Benchmarking/Unsupervised/Retina/output/stab_s4.csv",
          row.names = FALSE)
