######################### Random 1 ################################
library(StabMap)
library(scran)
library(Seurat)
library(patchwork)

setwd("./PBMC_TEA")
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

atac.list <- atac.list[c(1,4)]
rna.list <- rna.list[c(3,4)]
adt.list <- adt.list[c(2,4)]

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

assay_list = list(Batch1=d1,Batch2=d5,Batch3=d3,Batch4=rbind(d2,d4,d6))
gc()
mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("Batch4"),
               maxFeatures = (nrow(d1)+nrow(d3)+nrow(d5)),
               ncomponentsReference = 20,
               ncomponentsSubset = 20)
mm <- data.frame(batch = c(rep("B1",6194),rep("B2",6381),rep("B3",6412),rep("B4",6530)),
                 celltype = c(as.character(meta.list[[1]][,2]),
                              as.character(meta.list[[2]][,2]),as.character(meta.list[[3]][,2]),
                              as.character(meta.list[[4]][,2]))
)
rownames(mm) <- c(rownames(meta.list[[1]]),
                  rownames(meta.list[[2]]),
                  rownames(meta.list[[3]]),
                  rownames(meta.list[[4]]))
emb <- as.matrix(stab)
emb <- emb[rownames(mm),]

write.csv(as.matrix(emb),
          "./Benchmarking/Unsupervised/TEA_s1/output/stab_v1.csv",
          row.names = FALSE)

######################### Random 2 ################################
library(StabMap)
library(scran)
library(Seurat)
library(patchwork)

setwd("./PBMC_TEA")
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

mm <- data.frame(batch = c(rep("B1",6194),rep("B2",6381),rep("B3",6412),rep("B4",6530)),
                 celltype = c(as.character(meta.list[[1]][,2]),
                              as.character(meta.list[[2]][,2]),as.character(meta.list[[3]][,2]),
                              as.character(meta.list[[4]][,2]))
)
rownames(mm) <- c(rownames(meta.list[[1]]),
                  rownames(meta.list[[2]]),
                  rownames(meta.list[[3]]),
                  rownames(meta.list[[4]]))
emb <- as.matrix(stab)
emb <- emb[rownames(mm),]

write.csv(as.matrix(emb),
          "./Benchmarking/Unsupervised/TEA_s1/output/stab.csv",
          row.names = FALSE)

######################### Random 3 ################################
library(StabMap)
library(scran)
library(Seurat)
library(patchwork)

setwd("./PBMC_TEA")
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

atac.list <- atac.list[c(3,4)]
rna.list <- rna.list[c(1,3)]
adt.list <- adt.list[c(2,3)]

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

assay_list = list(Batch1=d3,
                  Batch2=d5,
                  Batch3=rbind(d1,d4,d6),
                  Batch4=d2)
#mf <- nrow(d1)
# rm(d1,d2,d3,d4,d5,d6)
gc()
mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("Batch3"),
               maxFeatures = (nrow(d1)+nrow(d3)+nrow(d5)),
               ncomponentsReference = 20,
               ncomponentsSubset = 20)

mm <- data.frame(batch = c(rep("B1",6194),rep("B2",6381),rep("B3",6412),rep("B4",6530)),
                 celltype = c(as.character(meta.list[[1]][,2]),
                              as.character(meta.list[[2]][,2]),as.character(meta.list[[3]][,2]),
                              as.character(meta.list[[4]][,2]))
)
rownames(mm) <- c(rownames(meta.list[[1]]),
                  rownames(meta.list[[2]]),
                  rownames(meta.list[[3]]),
                  rownames(meta.list[[4]]))
emb <- as.matrix(stab)
emb <- emb[rownames(mm),]

write.csv(as.matrix(emb),
          "./Benchmarking/Unsupervised/TEA_s1/output/stab_v2.csv",
          row.names = FALSE)

######################### Random 4 ################################
library(StabMap)
library(scran)
library(Seurat)
library(patchwork)

setwd("./PBMC_TEA")
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

atac.list <- atac.list[c(1,3)]
rna.list <- rna.list[c(1,4)]
adt.list <- adt.list[c(1,2)]

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

assay_list = list(Batch1=rbind(d1,d3,d5),
                  Batch2=d6,
                  Batch3=d2,
                  Batch4=d4)
mf <- nrow(d1)+nrow(d3)+nrow(d5)
rm(d1,d2,d3,d4,d5,d6)
gc()
mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("Batch1"),
               maxFeatures = mf,
               ncomponentsReference = 20,
               ncomponentsSubset = 20)

mm <- data.frame(batch = c(rep("B1",6194),rep("B2",6381),rep("B3",6412),rep("B4",6530)),
                 celltype = c(as.character(meta.list[[1]][,2]),
                              as.character(meta.list[[2]][,2]),as.character(meta.list[[3]][,2]),
                              as.character(meta.list[[4]][,2]))
)
rownames(mm) <- c(rownames(meta.list[[1]]),
                  rownames(meta.list[[2]]),
                  rownames(meta.list[[3]]),
                  rownames(meta.list[[4]]))
emb <- as.matrix(stab)
emb <- emb[rownames(mm),]

write.csv(as.matrix(emb),
          "./Benchmarking/Unsupervised/TEA_s1/output/stab_v3.csv",
          row.names = FALSE)

######################### Random 5 ################################
library(StabMap)
library(scran)
library(Seurat)
library(patchwork)

setwd("./PBMC_TEA")
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

atac.list <- atac.list[c(1,2)]
rna.list <- rna.list[c(1,4)]
adt.list <- adt.list[c(1,3)]

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

assay_list = list(Batch1=rbind(d1,d3,d5),
                  Batch2=d2,
                  Batch3=d6,
                  Batch4=d4)
mf <- nrow(d1)+nrow(d3)+nrow(d5)
rm(d1,d2,d3,d4,d5,d6)
gc()
mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("Batch1"),
               maxFeatures = mf,
               ncomponentsReference = 20,
               ncomponentsSubset = 20)

mm <- data.frame(batch = c(rep("B1",6194),rep("B2",6381),rep("B3",6412),rep("B4",6530)),
                 celltype = c(as.character(meta.list[[1]][,2]),
                              as.character(meta.list[[2]][,2]),as.character(meta.list[[3]][,2]),
                              as.character(meta.list[[4]][,2]))
)
rownames(mm) <- c(rownames(meta.list[[1]]),
                  rownames(meta.list[[2]]),
                  rownames(meta.list[[3]]),
                  rownames(meta.list[[4]]))
emb <- as.matrix(stab)
emb <- emb[rownames(mm),]

write.csv(as.matrix(emb),
          "./Benchmarking/Unsupervised/TEA_s1/output/stab_v4.csv",
          row.names = FALSE)