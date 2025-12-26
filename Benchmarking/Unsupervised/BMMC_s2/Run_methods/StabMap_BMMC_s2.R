######################### Random 1 ################################
library(StabMap)
library(scran)
library(Seurat)
library(patchwork)

setwd('./BMMC_CITE_Multiome')
path <- c('CITE/s1d1_cite','CITE/s1d2_cite','CITE/s1d3_cite',
          'Multiome/s1d1_multi','Multiome/s1d2_multi','Multiome/s1d3_multi')
rna.list <- list()
adt.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:6){
  
  meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
  rna.list[[path[i]]] <- readRDS(paste0(path[i],'/rna.rds'))
  
  if(i %in% c(1,2,3)){
    adt.list[[path[i]]] <- readRDS(paste0(path[i],'/adt.rds'))
  }
  if(i %in% c(4,5,6)){
    atac.list[[path[i]]] <- readRDS(paste0(path[i],'/atac.rds'))
  }
}

names(rna.list) <- NULL
names(adt.list) <- NULL
names(atac.list) <- NULL

for(i in 1:6){
  if(i %in% c(1,2,3)){
    adt.list[[i]] <- NormalizeData(adt.list[[i]])
    atac.list[[i]] <- NormalizeData(atac.list[[i]])
  }
  rna.list[[i]] <- NormalizeData(rna.list[[i]])
}
colnames(adt.list[[1]]) <- paste0(colnames(adt.list[[1]]),'-adt')
colnames(adt.list[[2]]) <- paste0(colnames(adt.list[[2]]),'-adt')

colnames(atac.list[[2]]) <- paste0(colnames(atac.list[[2]]),'-atac')
colnames(atac.list[[3]]) <- paste0(colnames(atac.list[[3]]),'-atac')

colnames(rna.list[[1]]) <- paste0(colnames(rna.list[[1]]),'-rna')
colnames(rna.list[[2]]) <- paste0(colnames(rna.list[[2]]),'-rna')
colnames(rna.list[[5]]) <- paste0(colnames(rna.list[[5]]),'-rna')
colnames(rna.list[[6]]) <- paste0(colnames(rna.list[[6]]),'-rna')



# rna for cite
rnac1 <- rna.list[[1]]
rnac2 <- rna.list[[2]]
rnac3 <- rna.list[[3]]

# rna for Multiome
rnam1 <- rna.list[[4]]
rnam2 <- rna.list[[5]]
rnam3 <- rna.list[[6]]

# adt for cite
adt1 <- adt.list[[1]]
adt2 <- adt.list[[2]]
adt3 <- adt.list[[3]]

# atac for Multiome
atac1 <- atac.list[[1]]
atac2 <- atac.list[[2]]
atac3 <- atac.list[[3]]


rm(rna.list,adt.list,atac.list)
gc()

rownames(rnac1) <- rownames(rnac2) <- rownames(rnac3) <- rownames(rnam1) <- rownames(rnam2) <- rownames(rnam3) <- paste0(rownames(rnac1),'-rna')
rownames(adt1) <- rownames(adt2) <- rownames(adt3) <- paste0(rownames(adt1),'-adt')
rownames(atac1) <- rownames(atac2) <- rownames(atac3) <- paste0(rownames(atac1),'-atac')

assay_list = list(CITE1rna=rnac1,
                  CITE1adt=adt1,
                  CITE2rna=rnac2,
                  CITE2adt=adt2,
                  CITE3=rbind(rnac3,adt3),
                  Multiome1=rbind(rnam1,atac1),
                  Multiome2rna=rnam2,
                  Multiome2atac=atac2,
                  Multiome3rna=rnam3,
                  Multiome3atac=atac3)
                  
rm(rnac1,rnac2,rnac3,rnam1,rnam2,rnam3,adt1,adt2,adt3,atac1,atac2,atac3)
gc()
mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("Multiome1"),
               maxFeatures = 1000000,
              ncomponentsReference = 20,
              ncomponentsSubset = 20)

meta <- readRDS('./Benchmarking/Unsupervised/BMMC_s2/output/meta_34.rds')

write.csv(as.matrix(stab[rownames(meta),]),
          "./Benchmarking/Unsupervised/BMMC_s2/output/stab_34.csv",
          row.names = FALSE)

######################### Random 2 ################################
library(StabMap)
library(scran)
library(Seurat)
library(patchwork)

setwd('./BMMC_CITE_Multiome')

path <- c('CITE/s1d1_cite','CITE/s1d2_cite','CITE/s1d3_cite',
          'Multiome/s1d1_multi','Multiome/s1d2_multi','Multiome/s1d3_multi')
rna.list <- list()
adt.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:6){
  
  meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
  rna.list[[path[i]]] <- readRDS(paste0(path[i],'/rna.rds'))
  
  if(i %in% c(1,2,3)){
    adt.list[[path[i]]] <- readRDS(paste0(path[i],'/adt.rds'))
  }
  if(i %in% c(4,5,6)){
    atac.list[[path[i]]] <- readRDS(paste0(path[i],'/atac.rds'))
  }
}

names(rna.list) <- NULL
names(adt.list) <- NULL
names(atac.list) <- NULL

for(i in 1:6){
  if(i %in% c(1,2,3)){
    adt.list[[i]] <- NormalizeData(adt.list[[i]])
    atac.list[[i]] <- NormalizeData(atac.list[[i]])
  }
  rna.list[[i]] <- NormalizeData(rna.list[[i]])
}

colnames(adt.list[[1]]) <- paste0(colnames(adt.list[[1]]),'-adt')
colnames(adt.list[[2]]) <- paste0(colnames(adt.list[[2]]),'-adt')

colnames(atac.list[[1]]) <- paste0(colnames(atac.list[[1]]),'-atac')
colnames(atac.list[[3]]) <- paste0(colnames(atac.list[[3]]),'-atac')

colnames(rna.list[[1]]) <- paste0(colnames(rna.list[[1]]),'-rna')
colnames(rna.list[[2]]) <- paste0(colnames(rna.list[[2]]),'-rna')
colnames(rna.list[[4]]) <- paste0(colnames(rna.list[[4]]),'-rna')
colnames(rna.list[[6]]) <- paste0(colnames(rna.list[[6]]),'-rna')

# rna for cite
rnac1 <- rna.list[[1]]
rnac2 <- rna.list[[2]]
rnac3 <- rna.list[[3]]

# rna for Multiome
rnam1 <- rna.list[[4]]
rnam2 <- rna.list[[5]]
rnam3 <- rna.list[[6]]

# adt for cite
adt1 <- adt.list[[1]]
adt2 <- adt.list[[2]]
adt3 <- adt.list[[3]]

# atac for Multiome
atac1 <- atac.list[[1]]
atac2 <- atac.list[[2]]
atac3 <- atac.list[[3]]


# rm(rna.list,adt.list,atac.list)
gc()

rownames(rnac1) <- rownames(rnac2) <- rownames(rnac3) <- rownames(rnam1) <- rownames(rnam2) <- rownames(rnam3) <- paste0(rownames(rnac1),'-rna')
rownames(adt1) <- rownames(adt2) <- rownames(adt3) <- paste0(rownames(adt1),'-adt')
rownames(atac1) <- rownames(atac2) <- rownames(atac3) <- paste0(rownames(atac1),'-atac')

assay_list = list(CITE1rna=rnac1,
                  CITE1adt=adt1,
                  CITE2rna=rnac2,
                  CITE2adt=adt2,
                  CITE3=rbind(rnac3,adt3),
                  Multiome1rna=rnam1,
                  Multiome1atac=atac1,
                  Multiome2=rbind(rnam2,atac2),
                  Multiome3rna=rnam3,
                  Multiome3atac=atac3)

rm(rnac1,rnac2,rnac3,rnam1,rnam2,rnam3,adt1,adt2,adt3,atac1,atac2,atac3)
gc()
mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("Multiome2"),
               maxFeatures = 1000000,
               ncomponentsReference = 20,
               ncomponentsSubset = 20)
meta_c1rna <- meta.list[[1]]
meta_c1rna$modality <- 'RNA'
meta_c1rna$batch <- 'CITE1_RNA'

meta_c1adt <- meta.list[[1]]
meta_c1adt$modality <- 'ADT'
meta_c1adt$batch <- 'CITE1_ADT'

meta_c2rna <- meta.list[[2]]
meta_c2rna$modality <- 'RNA'
meta_c2rna$batch <- 'CITE2_RNA'

meta_c2adt <- meta.list[[2]]
meta_c2adt$modality <- 'ADT'
meta_c2adt$batch <- 'CITE2_ADT'

meta_c3 <- meta.list[[3]]
meta_c3$modality <- 'CITE'
meta_c3$batch <- 'CITE3'

meta_m1rna <- meta.list[[4]]
meta_m1rna$modality <- 'RNA'
meta_m1rna$batch <- 'Multiome1_RNA'

meta_m1atac <- meta.list[[4]]
meta_m1atac$modality <- 'ATAC'
meta_m1atac$batch <- 'Multiome1_ATAC'

meta_m2 <- meta.list[[5]]
meta_m2$modality <- 'Multiome'
meta_m2$batch <- 'Multiome2'

meta_m3rna <- meta.list[[6]]
meta_m3rna$modality <- 'RNA'
meta_m3rna$batch <- 'Multiome3_RNA'

meta_m3atac <- meta.list[[6]]
meta_m3atac$modality <- 'ATAC'
meta_m3atac$batch <- 'Multiome3_ATAC'

rownames(meta_c1rna) <- colnames(rna.list[[1]])
rownames(meta_c1adt) <- colnames(adt.list[[1]])
rownames(meta_c2rna) <- colnames(rna.list[[2]])
rownames(meta_c2adt) <- colnames(adt.list[[2]])
rownames(meta_c3) <- colnames(rna.list[[3]])

rownames(meta_m1rna) <- colnames(rna.list[[4]])
rownames(meta_m1atac) <- colnames(atac.list[[1]])
rownames(meta_m2) <- colnames(rna.list[[5]])
rownames(meta_m3rna) <- colnames(rna.list[[6]])
rownames(meta_m3atac) <- colnames(atac.list[[3]])

meta <- rbind(meta_c1rna,meta_c1adt,
              meta_c2rna,meta_c2adt,
              meta_c3,
              meta_m1rna,meta_m1atac,
              meta_m2,
              meta_m3rna,meta_m3atac)

emb <- as.matrix(stab)
emb <- emb[rownames(meta),]

write.csv(emb,
          "./Benchmarking/Unsupervised/BMMC_s2/output/stab_35.csv",
          row.names = FALSE)

######################### Random 3 ################################
library(StabMap)
library(scran)
library(Seurat)
library(patchwork)

setwd('./BMMC_CITE_Multiome')

path <- c('CITE/s1d1_cite','CITE/s1d2_cite','CITE/s1d3_cite',
          'Multiome/s1d1_multi','Multiome/s1d2_multi','Multiome/s1d3_multi')
rna.list <- list()
adt.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:6){
  
  meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
  rna.list[[path[i]]] <- readRDS(paste0(path[i],'/rna.rds'))
  
  if(i %in% c(1,2,3)){
    adt.list[[path[i]]] <- readRDS(paste0(path[i],'/adt.rds'))
  }
  if(i %in% c(4,5,6)){
    atac.list[[path[i]]] <- readRDS(paste0(path[i],'/atac.rds'))
  }
}

names(rna.list) <- NULL
names(adt.list) <- NULL
names(atac.list) <- NULL

for(i in 1:6){
  if(i %in% c(1,2,3)){
    adt.list[[i]] <- NormalizeData(adt.list[[i]])
    atac.list[[i]] <- NormalizeData(atac.list[[i]])
  }
  rna.list[[i]] <- NormalizeData(rna.list[[i]])
}

colnames(adt.list[[1]]) <- paste0(colnames(adt.list[[1]]),'-adt')
colnames(adt.list[[2]]) <- paste0(colnames(adt.list[[2]]),'-adt')

colnames(atac.list[[1]]) <- paste0(colnames(atac.list[[1]]),'-atac')
colnames(atac.list[[2]]) <- paste0(colnames(atac.list[[2]]),'-atac')

colnames(rna.list[[1]]) <- paste0(colnames(rna.list[[1]]),'-rna')
colnames(rna.list[[2]]) <- paste0(colnames(rna.list[[2]]),'-rna')
colnames(rna.list[[4]]) <- paste0(colnames(rna.list[[4]]),'-rna')
colnames(rna.list[[5]]) <- paste0(colnames(rna.list[[5]]),'-rna')

# rna for cite
rnac1 <- rna.list[[1]]
rnac2 <- rna.list[[2]]
rnac3 <- rna.list[[3]]

# rna for Multiome
rnam1 <- rna.list[[4]]
rnam2 <- rna.list[[5]]
rnam3 <- rna.list[[6]]

# adt for cite
adt1 <- adt.list[[1]]
adt2 <- adt.list[[2]]
adt3 <- adt.list[[3]]

# atac for Multiome
atac1 <- atac.list[[1]]
atac2 <- atac.list[[2]]
atac3 <- atac.list[[3]]


# rm(rna.list,adt.list,atac.list)
gc()

rownames(rnac1) <- rownames(rnac2) <- rownames(rnac3) <- rownames(rnam1) <- rownames(rnam2) <- rownames(rnam3) <- paste0(rownames(rnac1),'-rna')
rownames(adt1) <- rownames(adt2) <- rownames(adt3) <- paste0(rownames(adt1),'-adt')
rownames(atac1) <- rownames(atac2) <- rownames(atac3) <- paste0(rownames(atac1),'-atac')

assay_list = list(CITE1rna=rnac1,
                  CITE1adt=adt1,
                  CITE2rna=rnac2,
                  CITE2adt=adt2,
                  CITE3=rbind(rnac3,adt3),
                  Multiome1rna=rnam1,
                  Multiome1atac=atac1,
                  Multiome2rna=rnam2,
                  Multiome2atac=atac2,
                  Multiome3=rbind(rnam3,atac3))

rm(rnac1,rnac2,rnac3,rnam1,rnam2,rnam3,adt1,adt2,adt3,atac1,atac2,atac3)
gc()
mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("Multiome3"),
               maxFeatures = 1000000,
               ncomponentsReference = 20,
               ncomponentsSubset = 20)

meta_c1rna <- meta.list[[1]]
meta_c1rna$modality <- 'RNA'
meta_c1rna$batch <- 'CITE1_RNA'

meta_c1adt <- meta.list[[1]]
meta_c1adt$modality <- 'ADT'
meta_c1adt$batch <- 'CITE1_ADT'

meta_c2rna <- meta.list[[2]]
meta_c2rna$modality <- 'RNA'
meta_c2rna$batch <- 'CITE2_RNA'

meta_c2adt <- meta.list[[2]]
meta_c2adt$modality <- 'ADT'
meta_c2adt$batch <- 'CITE2_ADT'

meta_c3 <- meta.list[[3]]
meta_c3$modality <- 'CITE'
meta_c3$batch <- 'CITE3'

meta_m1rna <- meta.list[[4]]
meta_m1rna$modality <- 'RNA'
meta_m1rna$batch <- 'Multiome1_RNA'

meta_m1atac <- meta.list[[4]]
meta_m1atac$modality <- 'ATAC'
meta_m1atac$batch <- 'Multiome1_ATAC'

meta_m2rna <- meta.list[[5]]
meta_m2rna$modality <- 'RNA'
meta_m2rna$batch <- 'Multiome2_RNA'

meta_m2atac <- meta.list[[5]]
meta_m2atac$modality <- 'ATAC'
meta_m2atac$batch <- 'Multiome2_ATAC'

meta_m3 <- meta.list[[6]]
meta_m3$modality <- 'Multiome'
meta_m3$batch <- 'Multiome3'

rownames(meta_c1rna) <- colnames(rna.list[[1]])
rownames(meta_c1adt) <- colnames(adt.list[[1]])
rownames(meta_c2rna) <- colnames(rna.list[[2]])
rownames(meta_c2adt) <- colnames(adt.list[[2]])
rownames(meta_c3) <- colnames(rna.list[[3]])

rownames(meta_m1rna) <- colnames(rna.list[[4]])
rownames(meta_m1atac) <- colnames(atac.list[[1]])
rownames(meta_m2rna) <- colnames(rna.list[[5]])
rownames(meta_m2atac) <- colnames(atac.list[[2]])
rownames(meta_m3) <- colnames(rna.list[[6]])

meta <- rbind(meta_c1rna,meta_c1adt,
              meta_c2rna,meta_c2adt,
              meta_c3,
              meta_m1rna,meta_m1atac,
              meta_m2rna,meta_m2atac,
              meta_m3)

emb <- as.matrix(stab)
emb <- emb[rownames(meta),]

write.csv(emb,
          "./Benchmarking/Unsupervised/BMMC_s2/output/stab_36.csv",
          row.names = FALSE)

######################### Random 4 ################################
library(StabMap)
library(scran)
library(Seurat)
library(patchwork)

setwd('./BMMC_CITE_Multiome/')

path <- c('CITE/s1d1_cite','CITE/s1d2_cite','CITE/s1d3_cite',
          'Multiome/s1d1_multi','Multiome/s1d2_multi','Multiome/s1d3_multi')
rna.list <- list()
adt.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:6){
  
  meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
  rna.list[[path[i]]] <- readRDS(paste0(path[i],'/rna.rds'))
  
  if(i %in% c(1,2,3)){
    adt.list[[path[i]]] <- readRDS(paste0(path[i],'/adt.rds'))
  }
  if(i %in% c(4,5,6)){
    atac.list[[path[i]]] <- readRDS(paste0(path[i],'/atac.rds'))
  }
}

names(rna.list) <- NULL
names(adt.list) <- NULL
names(atac.list) <- NULL

for(i in 1:6){
  if(i %in% c(1,2,3)){
    adt.list[[i]] <- NormalizeData(adt.list[[i]])
    atac.list[[i]] <- NormalizeData(atac.list[[i]])
  }
  rna.list[[i]] <- NormalizeData(rna.list[[i]])
}

colnames(adt.list[[1]]) <- paste0(colnames(adt.list[[1]]),'-adt')
colnames(adt.list[[3]]) <- paste0(colnames(adt.list[[3]]),'-adt')

colnames(atac.list[[2]]) <- paste0(colnames(atac.list[[2]]),'-atac')
colnames(atac.list[[3]]) <- paste0(colnames(atac.list[[3]]),'-atac')

colnames(rna.list[[1]]) <- paste0(colnames(rna.list[[1]]),'-rna')
colnames(rna.list[[3]]) <- paste0(colnames(rna.list[[3]]),'-rna')
colnames(rna.list[[5]]) <- paste0(colnames(rna.list[[5]]),'-rna')
colnames(rna.list[[6]]) <- paste0(colnames(rna.list[[6]]),'-rna')

# rna for cite
rnac1 <- rna.list[[1]]
rnac2 <- rna.list[[2]]
rnac3 <- rna.list[[3]]

# rna for Multiome
rnam1 <- rna.list[[4]]
rnam2 <- rna.list[[5]]
rnam3 <- rna.list[[6]]

# adt for cite
adt1 <- adt.list[[1]]
adt2 <- adt.list[[2]]
adt3 <- adt.list[[3]]

# atac for Multiome
atac1 <- atac.list[[1]]
atac2 <- atac.list[[2]]
atac3 <- atac.list[[3]]


# rm(rna.list,adt.list,atac.list)
gc()

rownames(rnac1) <- rownames(rnac2) <- rownames(rnac3) <- rownames(rnam1) <- rownames(rnam2) <- rownames(rnam3) <- paste0(rownames(rnac1),'-rna')
rownames(adt1) <- rownames(adt2) <- rownames(adt3) <- paste0(rownames(adt1),'-adt')
rownames(atac1) <- rownames(atac2) <- rownames(atac3) <- paste0(rownames(atac1),'-atac')

assay_list = list(CITE1rna=rnac1,
                  CITE1adt=adt1,
                  CITE2=rbind(rnac2,adt2),
                  CITE3rna=rnac3,
                  CITE3adt=adt3,
                  Multiome1=rbind(rnam1,atac1),
                  Multiome2rna=rnam2,
                  Multiome2atac=atac2,
                  Multiome3rna=rnam3,
                  Multiome3atac=atac3)

rm(rnac1,rnac2,rnac3,rnam1,rnam2,rnam3,adt1,adt2,adt3,atac1,atac2,atac3)
gc()
mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("Multiome1"),
               maxFeatures = 1000000,
               ncomponentsReference = 20,
               ncomponentsSubset = 20)

meta_c1rna <- meta.list[[1]]
meta_c1rna$modality <- 'RNA'
meta_c1rna$batch <- 'CITE1_RNA'

meta_c1adt <- meta.list[[1]]
meta_c1adt$modality <- 'ADT'
meta_c1adt$batch <- 'CITE1_ADT'

meta_c2 <- meta.list[[2]]
meta_c2$modality <- 'CITE'
meta_c2$batch <- 'CITE2'

meta_c3rna <- meta.list[[3]]
meta_c3rna$modality <- 'RNA'
meta_c3rna$batch <- 'CITE3_RNA'

meta_c3adt <- meta.list[[3]]
meta_c3adt$modality <- 'ADT'
meta_c3adt$batch <- 'CITE3_ADT'


meta_m1 <- meta.list[[4]]
meta_m1$modality <- 'Multiome'
meta_m1$batch <- 'Multiome1'


meta_m2rna <- meta.list[[5]]
meta_m2rna$modality <- 'RNA'
meta_m2rna$batch <- 'Multiome2_RNA'

meta_m2atac <- meta.list[[5]]
meta_m2atac$modality <- 'ATAC'
meta_m2atac$batch <- 'Multiome2_ATAC'

meta_m3rna <- meta.list[[6]]
meta_m3rna$modality <- 'RNA'
meta_m3rna$batch <- 'Multiome3_RNA'

meta_m3atac <- meta.list[[6]]
meta_m3atac$modality <- 'ATAC'
meta_m3atac$batch <- 'Multiome3_ATAC'

rownames(meta_c1rna) <- colnames(rna.list[[1]])
rownames(meta_c1adt) <- colnames(adt.list[[1]])
rownames(meta_c2) <- colnames(rna.list[[2]])
rownames(meta_c3rna) <- colnames(rna.list[[3]])
rownames(meta_c3adt) <- colnames(adt.list[[3]])

rownames(meta_m1) <- colnames(rna.list[[4]])
rownames(meta_m2atac) <- colnames(atac.list[[2]])
rownames(meta_m2rna) <- colnames(rna.list[[5]])
rownames(meta_m3atac) <- colnames(atac.list[[3]])
rownames(meta_m3rna) <- colnames(rna.list[[6]])

meta <- rbind(meta_c1rna,meta_c1adt,
              meta_c2,
              meta_c3rna,meta_c3adt,
              meta_m1,
              meta_m2rna,meta_m2atac,
              meta_m3rna,meta_m3atac)
emb <- as.matrix(stab)
emb <- emb[rownames(meta),]

stab_umap1 <- uwot::umap(as.matrix(emb),min_dist = 0.3)
dim.plot0(stab_umap1,meta,"batch")
dim.plot0(stab_umap1,meta,"celltype.l1")
dim.plot0(stab_umap1,meta,"celltype.l2")

write.csv(emb,
          "./Benchmarking/Unsupervised/BMMC_s2/output/stab_24.csv",
          row.names = FALSE)


######################### Random 5 ################################
library(StabMap)
library(scran)
library(Seurat)
library(patchwork)

setwd('./BMMC_CITE_Multiome')

path <- c('CITE/s1d1_cite','CITE/s1d2_cite','CITE/s1d3_cite',
          'Multiome/s1d1_multi','Multiome/s1d2_multi','Multiome/s1d3_multi')
rna.list <- list()
adt.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:6){
  
  meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
  rna.list[[path[i]]] <- readRDS(paste0(path[i],'/rna.rds'))
  
  if(i %in% c(1,2,3)){
    adt.list[[path[i]]] <- readRDS(paste0(path[i],'/adt.rds'))
  }
  if(i %in% c(4,5,6)){
    atac.list[[path[i]]] <- readRDS(paste0(path[i],'/atac.rds'))
  }
}

names(rna.list) <- NULL
names(adt.list) <- NULL
names(atac.list) <- NULL

for(i in 1:6){
  if(i %in% c(1,2,3)){
    adt.list[[i]] <- NormalizeData(adt.list[[i]])
    atac.list[[i]] <- NormalizeData(atac.list[[i]])
  }
  rna.list[[i]] <- NormalizeData(rna.list[[i]])
}

colnames(adt.list[[1]]) <- paste0(colnames(adt.list[[1]]),'-adt')
colnames(adt.list[[3]]) <- paste0(colnames(adt.list[[3]]),'-adt')

colnames(atac.list[[1]]) <- paste0(colnames(atac.list[[1]]),'-atac')
colnames(atac.list[[2]]) <- paste0(colnames(atac.list[[2]]),'-atac')

colnames(rna.list[[1]]) <- paste0(colnames(rna.list[[1]]),'-rna')
colnames(rna.list[[3]]) <- paste0(colnames(rna.list[[3]]),'-rna')
colnames(rna.list[[4]]) <- paste0(colnames(rna.list[[4]]),'-rna')
colnames(rna.list[[5]]) <- paste0(colnames(rna.list[[5]]),'-rna')

# rna for cite
rnac1 <- rna.list[[1]]
rnac2 <- rna.list[[2]]
rnac3 <- rna.list[[3]]

# rna for Multiome
rnam1 <- rna.list[[4]]
rnam2 <- rna.list[[5]]
rnam3 <- rna.list[[6]]

# adt for cite
adt1 <- adt.list[[1]]
adt2 <- adt.list[[2]]
adt3 <- adt.list[[3]]

# atac for Multiome
atac1 <- atac.list[[1]]
atac2 <- atac.list[[2]]
atac3 <- atac.list[[3]]


# rm(rna.list,adt.list,atac.list)
gc()

rownames(rnac1) <- rownames(rnac2) <- rownames(rnac3) <- rownames(rnam1) <- rownames(rnam2) <- rownames(rnam3) <- paste0(rownames(rnac1),'-rna')
rownames(adt1) <- rownames(adt2) <- rownames(adt3) <- paste0(rownames(adt1),'-adt')
rownames(atac1) <- rownames(atac2) <- rownames(atac3) <- paste0(rownames(atac1),'-atac')

assay_list = list(CITE1rna=rnac1,
                  CITE1adt=adt1,
                  CITE2=rbind(rnac2,adt2),
                  CITE3rna=rnac3,
                  CITE3adt=adt3,
                  Multiome1rna=rnam1,
                  Multiome1atac=atac1,
                  Multiome2rna=rnam2,
                  Multiome2atac=atac2,
                  Multiome3=rbind(rnam3,atac3))

rm(rnac1,rnac2,rnac3,rnam1,rnam2,rnam3,adt1,adt2,adt3,atac1,atac2,atac3)
gc()
mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("Multiome3"),
               maxFeatures = 1000000,
               ncomponentsReference = 20,
               ncomponentsSubset = 20)

meta_c1rna <- meta.list[[1]]
meta_c1rna$modality <- 'RNA'
meta_c1rna$batch <- 'CITE1_RNA'

meta_c1adt <- meta.list[[1]]
meta_c1adt$modality <- 'ADT'
meta_c1adt$batch <- 'CITE1_ADT'

meta_c2 <- meta.list[[2]]
meta_c2$modality <- 'CITE'
meta_c2$batch <- 'CITE2'

meta_c3rna <- meta.list[[3]]
meta_c3rna$modality <- 'RNA'
meta_c3rna$batch <- 'CITE3_RNA'

meta_c3adt <- meta.list[[3]]
meta_c3adt$modality <- 'ADT'
meta_c3adt$batch <- 'CITE3_ADT'


meta_m1rna <- meta.list[[4]]
meta_m1rna$modality <- 'RNA'
meta_m1rna$batch <- 'Multiome1_RNA'

meta_m1atac <- meta.list[[4]]
meta_m1atac$modality <- 'ATAC'
meta_m1atac$batch <- 'Multiome1_ATAC'

meta_m2rna <- meta.list[[5]]
meta_m2rna$modality <- 'RNA'
meta_m2rna$batch <- 'Multiome2_RNA'

meta_m2atac <- meta.list[[5]]
meta_m2atac$modality <- 'ATAC'
meta_m2atac$batch <- 'Multiome2_ATAC'

meta_m3 <- meta.list[[6]]
meta_m3$modality <- 'Multiome'
meta_m3$batch <- 'Multiome3'

rownames(meta_c1rna) <- colnames(rna.list[[1]])
rownames(meta_c1adt) <- colnames(adt.list[[1]])
rownames(meta_c2) <- colnames(rna.list[[2]])
rownames(meta_c3rna) <- colnames(rna.list[[3]])
rownames(meta_c3adt) <- colnames(adt.list[[3]])

rownames(meta_m1rna) <- colnames(rna.list[[4]])
rownames(meta_m1atac) <- colnames(atac.list[[1]])
rownames(meta_m2rna) <- colnames(rna.list[[5]])
rownames(meta_m2atac) <- colnames(atac.list[[2]])
rownames(meta_m3) <- colnames(rna.list[[6]])

meta <- rbind(meta_c1rna,meta_c1adt,
              meta_c2,
              meta_c3rna,meta_c3adt,
              meta_m1rna,meta_m1atac,
              meta_m2rna,meta_m2atac,
              meta_m3)
emb <- as.matrix(stab)
emb <- emb[rownames(meta),]

stab_umap1 <- uwot::umap(as.matrix(emb),min_dist = 0.3)
dim.plot0(stab_umap1,meta,"batch")
dim.plot0(stab_umap1,meta,"celltype.l1")
dim.plot0(stab_umap1,meta,"celltype.l2")

write.csv(emb,
          "./Benchmarking/Unsupervised/BMMC_s2/output/stab_26.csv",
          row.names = FALSE)
