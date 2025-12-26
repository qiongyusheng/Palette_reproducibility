######################### Random 1 ################################
library(StabMap)
library(scran)
library(Seurat)
library(patchwork)
setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s1d1_cite','./CITE/s1d2_cite','./CITE/s1d3_cite')
rna1.list <- list()
adt.list <- list()
meta1.list <- list()
for(i in 1:3){
    rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
    adt.list[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))
    meta1.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}

path2 <- c('./Multiome/s1d1_multi','./Multiome/s1d3_multi')
rna2.list <- list()
atac.list <- list()
meta2.list <- list()
for(i in 1:2){
    rna2.list[[i]] <- readRDS(paste0(path2[i],'/rna.rds'))
    atac.list[[i]] <- readRDS(paste0(path2[i],'/atac.rds'))
    meta2.list[[i]] <- readRDS(paste0(path2[i],'/meta.rds'))
}

for(i in 1:3){
    rna1.list[[i]] <- NormalizeData(rna1.list[[i]])
    adt.list[[i]] <- NormalizeData(adt.list[[i]])
}
for(i in 1:2){
    rna2.list[[i]] <- NormalizeData(rna2.list[[i]])
    atac.list[[i]] <- NormalizeData(atac.list[[i]])
}


# rna for cite
rnac1 <- rna1.list[[1]]
rnac2 <- rna1.list[[2]]
rnac3 <- rna1.list[[3]]

# rna for Multiome
rnam1 <- rna2.list[[1]]
rnam2 <- rna2.list[[2]]

# adt for cite
adt1 <- adt.list[[1]]
adt2 <- adt.list[[2]]
adt3 <- adt.list[[3]]

# atac for Multiome
atac1 <- atac.list[[1]]
atac2 <- atac.list[[2]]

rm(rna1.list,rna2.list,adt.list,atac.list)
gc()

rownames(rnac1) <- rownames(rnac2) <- rownames(rnac3) <- rownames(rnam1) <- rownames(rnam2) <- paste0(rownames(rnac1),'-rna')
rownames(adt1) <- rownames(adt2) <- rownames(adt3) <- paste0(rownames(adt1),'-adt')
rownames(atac1) <- rownames(atac2) <- paste0(rownames(atac1),'-atac')

assay_list = list(CITE1=rbind(rnac1,adt1),
                  CITE2=rbind(rnac2,adt2),
                  CITE3=rbind(rnac3,adt3),
                  Multiome1=rbind(rnam1,atac1),
                  Multiome3=rbind(rnam2,atac2))

rm(rnac1,adt1,rnac2,adt2,rnac3,adt3,rnam1,atac1,rnam2,atac2)

mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("Multiome1"),
               maxFeatures = 1000000,
              ncomponentsReference = 20,
              ncomponentsSubset = 20)

mm <- rbind(meta1.list[[1]],meta1.list[[2]],meta1.list[[3]],
            meta2.list[[1]],meta2.list[[2]])
mm$batch <- c(rep("CITE1",5227),rep("CITE2",4978),rep("CITE3",6106),
              rep("Multiome1",6224),rep("Multiome3",4279))
emb <- stab[rownames(mm),]

write.csv(as.matrix(emb),
          "./Benchmarking/Unsupervised/BMMC_s1/output/stab_drop5.csv",
          row.names = FALSE)


######################### Random 2 ################################
library(StabMap)
library(scran)
library(Seurat)
library(patchwork)

setwd('./BMMC_CITE_Multiome')

path1 <- c('./CITE/s1d2_cite','./CITE/s1d3_cite')
rna1.list <- list()
adt.list <- list()
meta1.list <- list()
for(i in 1:2){
  rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  adt.list[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))
  meta1.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}

path2 <- c('./Multiome/s1d1_multi','./Multiome/s1d2_multi','./Multiome/s1d3_multi')
rna2.list <- list()
atac.list <- list()
meta2.list <- list()
for(i in 1:3){
  rna2.list[[i]] <- readRDS(paste0(path2[i],'/rna.rds'))
  atac.list[[i]] <- readRDS(paste0(path2[i],'/atac.rds'))
  meta2.list[[i]] <- readRDS(paste0(path2[i],'/meta.rds'))
}

for(i in 1:2){
  rna1.list[[i]] <- NormalizeData(rna1.list[[i]])
  adt.list[[i]] <- NormalizeData(adt.list[[i]])
}
for(i in 1:3){
  rna2.list[[i]] <- NormalizeData(rna2.list[[i]])
  atac.list[[i]] <- NormalizeData(atac.list[[i]])
}


# rna for cite
rnac1 <- rna1.list[[1]]
rnac2 <- rna1.list[[2]]

# rna for Multiome
rnam1 <- rna2.list[[1]]
rnam2 <- rna2.list[[2]]
rnam3 <- rna2.list[[3]]

# adt for cite
adt1 <- adt.list[[1]]
adt2 <- adt.list[[2]]

# atac for Multiome
atac1 <- atac.list[[1]]
atac2 <- atac.list[[2]]
atac3 <- atac.list[[3]]

rm(rna1.list,rna2.list,adt.list,atac.list)
gc()

rownames(rnac1) <- rownames(rnac2) <- rownames(rnam1) <- rownames(rnam2) <- rownames(rnam3) <- paste0(rownames(rnac1),'-rna')
rownames(adt1) <- rownames(adt2) <- paste0(rownames(adt1),'-adt')
rownames(atac1) <- rownames(atac2) <- rownames(atac3) <- paste0(rownames(atac1),'-atac')

assay_list = list(CITE2=rbind(rnac1,adt1),
                  CITE3=rbind(rnac2,adt2),
                  Multiome1=rbind(rnam1,atac1),
                  Multiome2=rbind(rnam2,atac2),
                  Multiome3=rbind(rnam3,atac3))

rm(rnac1,adt1,rnac2,adt2,rnam1,atac1,rnam2,atac2,rnam3,atac3)

mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("Multiome2"),
               maxFeatures = 1000000,
               ncomponentsReference = 20,
               ncomponentsSubset = 20)

mm <- rbind(meta1.list[[1]],meta1.list[[2]],
            meta2.list[[1]],meta2.list[[2]],meta2.list[[3]])
mm$batch <- c(rep("CITE2",4978),rep("CITE3",6106),
              rep("Multiome1",6224),rep("Multiome2",6740),rep("Multiome3",4279))
emb <- stab[rownames(mm),]

write.csv(as.matrix(emb),
          "./Benchmarking/Unsupervised/BMMC_s1/output/stab_drop1.csv",
          row.names = FALSE)


######################### Random 3 ################################
library(StabMap)
library(scran)
library(Seurat)
library(patchwork)

setwd('./BMMC_CITE_Multiome')

path1 <- c('./CITE/s1d1_cite','./CITE/s1d2_cite','./CITE/s1d3_cite')
rna1.list <- list()
adt.list <- list()
meta1.list <- list()
for(i in 1:3){
  rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  adt.list[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))
  meta1.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}

path2 <- c('./Multiome/s1d1_multi','./Multiome/s1d2_multi')
rna2.list <- list()
atac.list <- list()
meta2.list <- list()
for(i in 1:2){
  rna2.list[[i]] <- readRDS(paste0(path2[i],'/rna.rds'))
  atac.list[[i]] <- readRDS(paste0(path2[i],'/atac.rds'))
  meta2.list[[i]] <- readRDS(paste0(path2[i],'/meta.rds'))
}

for(i in 1:3){
  rna1.list[[i]] <- NormalizeData(rna1.list[[i]])
  adt.list[[i]] <- NormalizeData(adt.list[[i]])
}
for(i in 1:2){
  rna2.list[[i]] <- NormalizeData(rna2.list[[i]])
  atac.list[[i]] <- NormalizeData(atac.list[[i]])
}


# rna for cite
rnac1 <- rna1.list[[1]]
rnac2 <- rna1.list[[2]]
rnac3 <- rna1.list[[3]]

# rna for Multiome
rnam1 <- rna2.list[[1]]
rnam2 <- rna2.list[[2]]

# adt for cite
adt1 <- adt.list[[1]]
adt2 <- adt.list[[2]]
adt3 <- adt.list[[3]]

# atac for Multiome
atac1 <- atac.list[[1]]
atac2 <- atac.list[[2]]

rm(rna1.list,rna2.list,adt.list,atac.list)
gc()

rownames(rnac1) <- rownames(rnac2) <- rownames(rnac3) <- rownames(rnam1) <- rownames(rnam2) <- paste0(rownames(rnac1),'-rna')
rownames(adt1) <- rownames(adt2) <- rownames(adt3) <- paste0(rownames(adt1),'-adt')
rownames(atac1) <- rownames(atac2) <- paste0(rownames(atac1),'-atac')

assay_list = list(CITE1=rbind(rnac1,adt1),
                  CITE2=rbind(rnac2,adt2),
                  CITE3=rbind(rnac3,adt3),
                  Multiome1=rbind(rnam1,atac1),
                  Multiome2=rbind(rnam2,atac2))

rm(rnac1,adt1,rnac2,adt2,rnac3,adt3,rnam1,atac1,rnam2,atac2)

mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("Multiome2"),
               maxFeatures = 1000000,
               ncomponentsReference = 20,
               ncomponentsSubset = 20)
mm <- rbind(meta1.list[[1]],meta1.list[[2]],meta1.list[[3]],
            meta2.list[[1]],meta2.list[[2]])
mm$batch <- c(rep("CITE1",5227),rep("CITE2",4978),rep("CITE3",6106),
              rep("Multiome1",6224),rep("Multiome2",6740))
emb <- stab[rownames(mm),]

write.csv(as.matrix(emb),
          "./Benchmarking/Unsupervised/BMMC_s1/output/stab_drop6.csv",
          row.names = FALSE)


######################### Random 4 ################################
library(StabMap)
library(scran)
library(Seurat)
library(patchwork)

setwd('./BMMC_CITE_Multiome')

path1 <- c('./CITE/s1d1_cite','./CITE/s1d2_cite','./CITE/s1d3_cite')
rna1.list <- list()
adt.list <- list()
meta1.list <- list()
for(i in 1:3){
  rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  adt.list[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))
  meta1.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}

path2 <- c('./Multiome/s1d2_multi','./Multiome/s1d3_multi')
rna2.list <- list()
atac.list <- list()
meta2.list <- list()
for(i in 1:2){
  rna2.list[[i]] <- readRDS(paste0(path2[i],'/rna.rds'))
  atac.list[[i]] <- readRDS(paste0(path2[i],'/atac.rds'))
  meta2.list[[i]] <- readRDS(paste0(path2[i],'/meta.rds'))
}

for(i in 1:3){
  rna1.list[[i]] <- NormalizeData(rna1.list[[i]])
  adt.list[[i]] <- NormalizeData(adt.list[[i]])
}
for(i in 1:2){
  rna2.list[[i]] <- NormalizeData(rna2.list[[i]])
  atac.list[[i]] <- NormalizeData(atac.list[[i]])
}


# rna for cite
rnac1 <- rna1.list[[1]]
rnac2 <- rna1.list[[2]]
rnac3 <- rna1.list[[3]]

# rna for Multiome
rnam1 <- rna2.list[[1]]
rnam2 <- rna2.list[[2]]

# adt for cite
adt1 <- adt.list[[1]]
adt2 <- adt.list[[2]]
adt3 <- adt.list[[3]]

# atac for Multiome
atac1 <- atac.list[[1]]
atac2 <- atac.list[[2]]

rm(rna1.list,rna2.list,adt.list,atac.list)
gc()

rownames(rnac1) <- rownames(rnac2) <- rownames(rnac3) <- rownames(rnam1) <- rownames(rnam2) <- paste0(rownames(rnac1),'-rna')
rownames(adt1) <- rownames(adt2) <- rownames(adt3) <- paste0(rownames(adt1),'-adt')
rownames(atac1) <- rownames(atac2) <- paste0(rownames(atac1),'-atac')

assay_list = list(CITE1=rbind(rnac1,adt1),
                  CITE2=rbind(rnac2,adt2),
                  CITE3=rbind(rnac3,adt3),
                  Multiome2=rbind(rnam1,atac1),
                  Multiome3=rbind(rnam2,atac2))

rm(rnac1,adt1,rnac2,adt2,rnac3,adt3,rnam1,atac1,rnam2,atac2)

mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("Multiome2"),
               maxFeatures = 1000000,
               ncomponentsReference = 20,
               ncomponentsSubset = 20)
mm <- rbind(meta1.list[[1]],meta1.list[[2]],meta1.list[[3]],
            meta2.list[[1]],meta2.list[[2]])
mm$batch <- c(rep("CITE1",5227),rep("CITE2",4978),rep("CITE3",6106),
              rep("Multiome2",6740),rep("Multiome3",4279))
emb <- stab[rownames(mm),]

write.csv(as.matrix(emb),
          "./Benchmarking/Unsupervised/BMMC_s1/output/stab_drop4.csv",
          row.names = FALSE)


######################### Random 5 ################################
library(StabMap)
library(scran)
library(Seurat)
library(patchwork)

setwd('./BMMC_CITE_Multiome')

path1 <- c('./CITE/s1d1_cite','./CITE/s1d3_cite')
rna1.list <- list()
adt.list <- list()
meta1.list <- list()
for(i in 1:2){
  rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  adt.list[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))
  meta1.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}

path2 <- c('./Multiome/s1d1_multi','./Multiome/s1d2_multi','./Multiome/s1d3_multi')
rna2.list <- list()
atac.list <- list()
meta2.list <- list()
for(i in 1:3){
  rna2.list[[i]] <- readRDS(paste0(path2[i],'/rna.rds'))
  atac.list[[i]] <- readRDS(paste0(path2[i],'/atac.rds'))
  meta2.list[[i]] <- readRDS(paste0(path2[i],'/meta.rds'))
}

for(i in 1:2){
  rna1.list[[i]] <- NormalizeData(rna1.list[[i]])
  adt.list[[i]] <- NormalizeData(adt.list[[i]])
}
for(i in 1:3){
  rna2.list[[i]] <- NormalizeData(rna2.list[[i]])
  atac.list[[i]] <- NormalizeData(atac.list[[i]])
}


# rna for cite
rnac1 <- rna1.list[[1]]
rnac2 <- rna1.list[[2]]

# rna for Multiome
rnam1 <- rna2.list[[1]]
rnam2 <- rna2.list[[2]]
rnam3 <- rna2.list[[3]]

# adt for cite
adt1 <- adt.list[[1]]
adt2 <- adt.list[[2]]

# atac for Multiome
atac1 <- atac.list[[1]]
atac2 <- atac.list[[2]]
atac3 <- atac.list[[3]]

rm(rna1.list,rna2.list,adt.list,atac.list)
gc()

rownames(rnac1) <- rownames(rnac2) <- rownames(rnam1) <- rownames(rnam2) <- rownames(rnam3) <- paste0(rownames(rnac1),'-rna')
rownames(adt1) <- rownames(adt2) <- paste0(rownames(adt1),'-adt')
rownames(atac1) <- rownames(atac2) <- rownames(atac3) <- paste0(rownames(atac1),'-atac')

assay_list = list(CITE1=rbind(rnac1,adt1),
                  CITE3=rbind(rnac2,adt2),
                  Multiome1=rbind(rnam1,atac1),
                  Multiome2=rbind(rnam2,atac2),
                  Multiome3=rbind(rnam3,atac3))

rm(rnac1,adt1,rnac2,adt2,rnam1,atac1,rnam2,atac2,rnam3,atac3)

mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("Multiome2"),
               maxFeatures = 1000000,
               ncomponentsReference = 20,
               ncomponentsSubset = 20)

mm <- rbind(meta1.list[[1]],meta1.list[[2]],
            meta2.list[[1]],meta2.list[[2]],meta2.list[[3]])
emb <- stab[rownames(mm),]

write.csv(as.matrix(emb),
          "./Benchmarking/Unsupervised/BMMC_s1/output/stab_drop2.csv",
          row.names = FALSE)