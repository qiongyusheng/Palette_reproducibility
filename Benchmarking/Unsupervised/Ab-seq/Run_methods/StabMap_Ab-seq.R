######################### Random 1 ################################
library(StabMap)
library(scran)
library(Seurat)
library(patchwork)

setwd('./Young_Age')

paths <- c('Age1','Age2','Age3','BM1','BM2','BM3')
adt.list <- list()
rna.list <- list()
meta.list <- list()
for(i in paths){
    adt.list[[i]] <- readRDS(paste0('./',i,'/adt.rds'))
    rna.list[[i]] <- readRDS(paste0('./',i,'/rna.rds'))
    meta.list[[i]] <- readRDS(paste0('./',i,'/meta.rds'))
}
names(adt.list) <- paste0(names(adt.list),'-adt')
names(rna.list) <- paste0(names(rna.list),'-rna')
adt.list <- adt.list[c(1,3,4,5)]
rna.list <- rna.list[c(1,2,5,6)]

for(i in 1:4){
    rna.list[[i]] <- NormalizeData(rna.list[[i]])
    adt.list[[i]] <- NormalizeData(adt.list[[i]])
}

# ADT
adt_a1 <- adt.list[[1]]
adt_a3 <- adt.list[[2]]
adt_b1 <- adt.list[[3]]
adt_b2 <- adt.list[[4]]
# RNA
rna_a1 <- rna.list[[1]]
rna_a2 <- rna.list[[2]]
rna_b2 <- rna.list[[3]]
rna_b3 <- rna.list[[4]]
rm(rna.list,adt.list)

rownames(adt_a1) <- paste0(rownames(adt_a1),'-adt')
rownames(adt_a3) <- paste0(rownames(adt_a3),'-adt')
rownames(adt_b2) <- paste0(rownames(adt_b2),'-adt')
rownames(adt_b1) <- paste0(rownames(adt_b1),'-adt')

rownames(rna_a2) <- paste0(rownames(rna_a2),'-rna')
rownames(rna_a1) <- paste0(rownames(rna_a1),'-rna')
rownames(rna_b3) <- paste0(rownames(rna_b3),'-rna')
rownames(rna_b2) <- paste0(rownames(rna_b2),'-rna')


assay_list = list(Age1=rbind(rna_a1,adt_a1),Age2=rna_a2,Age3=adt_a3,
                  BM1=adt_b1,BM2=rbind(rna_b2,adt_b2),BM3=rna_b3)
mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("BM2"),
               maxFeatures = 10000,
              ncomponentsReference = 20,
              ncomponentsSubset = 20)

dim(stab)

cell_name <- unlist(lapply(meta.list,rownames))
stab <- stab[cell_name,]

out_dir <- './Benchmarking/Unsupervised/Ab-seq/output'
write.csv(as.matrix(stab),
          paste0(out_dir,'/stab_s1.csv'),
          row.names = FALSE)

######################### Random 2 ################################
library(StabMap)
library(scran)
library(Seurat)
library(patchwork)

setwd('./Young_Age')

paths <- c('Age1','Age2','Age3','BM1','BM2','BM3')
adt.list <- list()
rna.list <- list()
meta.list <- list()
for(i in paths){
  adt.list[[i]] <- readRDS(paste0('./',i,'/adt.rds'))
  rna.list[[i]] <- readRDS(paste0('./',i,'/rna.rds'))
  meta.list[[i]] <- readRDS(paste0('./',i,'/meta.rds'))
}
names(adt.list) <- paste0(names(adt.list),'-adt')
names(rna.list) <- paste0(names(rna.list),'-rna')
adt.list <- adt.list[c(1,3,5,6)]
rna.list <- rna.list[c(2,3,4,6)]

for(i in 1:4){
  rna.list[[i]] <- NormalizeData(rna.list[[i]])
  adt.list[[i]] <- NormalizeData(adt.list[[i]])
}

# ADT
adt_a1 <- adt.list[[1]]
adt_a3 <- adt.list[[2]]
adt_b2 <- adt.list[[3]]
adt_b3 <- adt.list[[4]]
# RNA
rna_a2 <- rna.list[[1]]
rna_a3 <- rna.list[[2]]
rna_b1 <- rna.list[[3]]
rna_b3 <- rna.list[[4]]
rm(rna.list,adt.list)

rownames(adt_a1) <- paste0(rownames(adt_a1),'-adt')
rownames(adt_a3) <- paste0(rownames(adt_a3),'-adt')
rownames(adt_b2) <- paste0(rownames(adt_b2),'-adt')
rownames(adt_b3) <- paste0(rownames(adt_b3),'-adt')

rownames(rna_a2) <- paste0(rownames(rna_a2),'-rna')
rownames(rna_a3) <- paste0(rownames(rna_a3),'-rna')
rownames(rna_b3) <- paste0(rownames(rna_b3),'-rna')
rownames(rna_b1) <- paste0(rownames(rna_b1),'-rna')

assay_list = list(Age1=adt_a1,Age2=rna_a2,Age3=rbind(rna_a3,adt_a3),
                  BM1=rna_b1,BM2=adt_b2,BM3=rbind(rna_b3,adt_b3))
mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("BM3"),
               maxFeatures = 10000,
               ncomponentsReference = 20,
               ncomponentsSubset = 20)

mm <- data.frame(batch = c(rep("Age1",6316),rep("Age2",4071),rep("Age3",10639),
                           rep("BM1",9751),rep("BM2",6963),rep("BM3",11317)),
                 celltype = c(as.character(meta.list[[1]][,4]),
                              as.character(meta.list[[2]][,4]),
                              as.character(meta.list[[3]][,4]),
                              as.character(meta.list[[4]][,4]),
                              as.character(meta.list[[5]][,4]),
                              as.character(meta.list[[6]][,4])),
                 celltype1 = c(as.character(meta.list[[1]][,5]),
                               as.character(meta.list[[2]][,5]),
                               as.character(meta.list[[3]][,5]),
                               as.character(meta.list[[4]][,5]),
                               as.character(meta.list[[5]][,5]),
                               as.character(meta.list[[6]][,5])),
                 celltype2 = c(as.character(meta.list[[1]][,3]),
                               as.character(meta.list[[2]][,3]),
                               as.character(meta.list[[3]][,3]),
                               as.character(meta.list[[4]][,3]),
                               as.character(meta.list[[5]][,3]),
                               as.character(meta.list[[6]][,3])),
                 modality = c(as.character(meta.list[[1]][,6]),
                              as.character(meta.list[[2]][,6]),
                              as.character(meta.list[[3]][,6]),
                              as.character(meta.list[[4]][,6]),
                              as.character(meta.list[[5]][,6]),
                              as.character(meta.list[[6]][,6]))
)

cell_name <- unlist(lapply(meta.list,rownames))
stab <- stab[cell_name,]

out_dir <- './Benchmarking/Unsupervised/Ab-seq/output'
write.csv(as.matrix(stab),
          paste0(out_dir,'/stab_s2.csv'),
          row.names = FALSE)


######################### Random 3 ################################
library(StabMap)
library(scran)
library(Seurat)
library(patchwork)

setwd('./Young_Age')

paths <- c('Age1','Age2','Age3','BM1','BM2','BM3')
adt.list <- list()
rna.list <- list()
meta.list <- list()
for(i in paths){
  adt.list[[i]] <- readRDS(paste0('./',i,'/adt.rds'))
  rna.list[[i]] <- readRDS(paste0('./',i,'/rna.rds'))
  meta.list[[i]] <- readRDS(paste0('./',i,'/meta.rds'))
}
names(adt.list) <- paste0(names(adt.list),'-adt')
names(rna.list) <- paste0(names(rna.list),'-rna')
adt.list <- adt.list[c(1,2,4,6)]
rna.list <- rna.list[c(2,3,4,5)]

for(i in 1:4){
  rna.list[[i]] <- NormalizeData(rna.list[[i]])
  adt.list[[i]] <- NormalizeData(adt.list[[i]])
}

# ADT
adt_a1 <- adt.list[[1]]
adt_a2 <- adt.list[[2]]
adt_b1 <- adt.list[[3]]
adt_b3 <- adt.list[[4]]
# RNA
rna_a2 <- rna.list[[1]]
rna_a3 <- rna.list[[2]]
rna_b1 <- rna.list[[3]]
rna_b2 <- rna.list[[4]]
rm(rna.list,adt.list)

rownames(adt_a1) <- paste0(rownames(adt_a1),'-adt')
rownames(adt_a2) <- paste0(rownames(adt_a2),'-adt')
rownames(adt_b1) <- paste0(rownames(adt_b1),'-adt')
rownames(adt_b3) <- paste0(rownames(adt_b3),'-adt')

rownames(rna_a2) <- paste0(rownames(rna_a2),'-rna')
rownames(rna_a3) <- paste0(rownames(rna_a3),'-rna')
rownames(rna_b1) <- paste0(rownames(rna_b1),'-rna')
rownames(rna_b2) <- paste0(rownames(rna_b2),'-rna')


assay_list = list(Age1=adt_a1,Age2=rbind(rna_a2,adt_a2),Age3=rna_a3,
                  BM1=rbind(rna_b1,adt_b1),BM2=rna_b2,BM3=adt_b3)
mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("BM1"),
               maxFeatures = 10000,
               ncomponentsReference = 20,
               ncomponentsSubset = 20)

mm <- data.frame(batch = c(rep("Age1",6316),rep("Age2",4071),rep("Age3",10639),
                           rep("BM1",9751),rep("BM2",6963),rep("BM3",11317)),
                 celltype = c(as.character(meta.list[[1]][,4]),
                              as.character(meta.list[[2]][,4]),
                              as.character(meta.list[[3]][,4]),
                              as.character(meta.list[[4]][,4]),
                              as.character(meta.list[[5]][,4]),
                              as.character(meta.list[[6]][,4])),
                 celltype1 = c(as.character(meta.list[[1]][,5]),
                               as.character(meta.list[[2]][,5]),
                               as.character(meta.list[[3]][,5]),
                               as.character(meta.list[[4]][,5]),
                               as.character(meta.list[[5]][,5]),
                               as.character(meta.list[[6]][,5])),
                 celltype2 = c(as.character(meta.list[[1]][,3]),
                               as.character(meta.list[[2]][,3]),
                               as.character(meta.list[[3]][,3]),
                               as.character(meta.list[[4]][,3]),
                               as.character(meta.list[[5]][,3]),
                               as.character(meta.list[[6]][,3])),
                 modality = c(as.character(meta.list[[1]][,6]),
                              as.character(meta.list[[2]][,6]),
                              as.character(meta.list[[3]][,6]),
                              as.character(meta.list[[4]][,6]),
                              as.character(meta.list[[5]][,6]),
                              as.character(meta.list[[6]][,6]))
)

cell_name <- unlist(lapply(meta.list,rownames))

stab <- stab[cell_name,]

out_dir <- './Benchmarking/Unsupervised/Ab-seq/output'
write.csv(as.matrix(stab),
          paste0(out_dir,'/stab_s3.csv'),
          row.names = FALSE)

######################### Random 4 ################################
library(StabMap)
library(scran)
library(Seurat)
library(patchwork)

setwd('./Young_Age')
out_dir <- './Benchmarking/Unsupervised/Ab-seq/output'

paths <- c('Age1','Age2','Age3','BM1','BM2','BM3')
adt.list <- list()
rna.list <- list()
meta.list <- list()
for(i in paths){
  adt.list[[i]] <- readRDS(paste0('./',i,'/adt.rds'))
  rna.list[[i]] <- readRDS(paste0('./',i,'/rna.rds'))
  meta.list[[i]] <- readRDS(paste0('./',i,'/meta.rds'))
}
names(adt.list) <- paste0(names(adt.list),'-adt')
names(rna.list) <- paste0(names(rna.list),'-rna')
adt.list <- adt.list[c(1,2,4,6)]
rna.list <- rna.list[c(1,3,4,5)]

for(i in 1:4){
  rna.list[[i]] <- NormalizeData(rna.list[[i]])
  adt.list[[i]] <- NormalizeData(adt.list[[i]])
}

# ADT
adt_a1 <- adt.list[[1]]
adt_a2 <- adt.list[[2]]
adt_b1 <- adt.list[[3]]
adt_b3 <- adt.list[[4]]
# RNA
rna_a1 <- rna.list[[1]]
rna_a3 <- rna.list[[2]]
rna_b1 <- rna.list[[3]]
rna_b2 <- rna.list[[4]]
rm(rna.list,adt.list)

rownames(adt_a1) <- paste0(rownames(adt_a1),'-adt')
rownames(adt_a2) <- paste0(rownames(adt_a2),'-adt')
rownames(adt_b1) <- paste0(rownames(adt_b1),'-adt')
rownames(adt_b3) <- paste0(rownames(adt_b3),'-adt')

rownames(rna_a1) <- paste0(rownames(rna_a1),'-rna')
rownames(rna_a3) <- paste0(rownames(rna_a3),'-rna')
rownames(rna_b1) <- paste0(rownames(rna_b1),'-rna')
rownames(rna_b2) <- paste0(rownames(rna_b2),'-rna')

assay_list = list(Age1=rbind(rna_a1,adt_a1),Age2=adt_a2,Age3=rna_a3,
                  BM1=rbind(rna_b1,adt_b1),BM2=rna_b2,BM3=adt_b3)
mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("BM1"),
               maxFeatures = 10000,
               ncomponentsReference = 50,
               ncomponentsSubset = 50)

meta <- Reduce(rbind,meta.list)
stab <- stab[rownames(stab),]

write.csv(as.matrix(stab),
          paste0(out_dir,'/stab.csv'),
          row.names = FALSE)

######################### Random 5 ################################

library(StabMap)
library(scran)
library(Seurat)
library(patchwork)

setwd('./Young_Age')

paths <- c('Age1','Age2','Age3','BM1','BM2','BM3')
adt.list <- list()
rna.list <- list()
meta.list <- list()
for(i in paths){
  adt.list[[i]] <- readRDS(paste0('./',i,'/adt.rds'))
  rna.list[[i]] <- readRDS(paste0('./',i,'/rna.rds'))
  meta.list[[i]] <- readRDS(paste0('./',i,'/meta.rds'))
}
names(adt.list) <- paste0(names(adt.list),'-adt')
names(rna.list) <- paste0(names(rna.list),'-rna')
adt.list <- adt.list[c(1,2,4,5)]
rna.list <- rna.list[c(1,3,5,6)]

for(i in 1:4){
  rna.list[[i]] <- NormalizeData(rna.list[[i]])
  adt.list[[i]] <- NormalizeData(adt.list[[i]])
}

# ADT
adt_a1 <- adt.list[[1]]
adt_a2 <- adt.list[[2]]
adt_b1 <- adt.list[[3]]
adt_b2 <- adt.list[[4]]
# RNA
rna_a1 <- rna.list[[1]]
rna_a3 <- rna.list[[2]]
rna_b2 <- rna.list[[3]]
rna_b3 <- rna.list[[4]]
rm(rna.list,adt.list)

rownames(adt_a1) <- paste0(rownames(adt_a1),'-adt')
rownames(adt_a2) <- paste0(rownames(adt_a2),'-adt')
rownames(adt_b1) <- paste0(rownames(adt_b1),'-adt')
rownames(adt_b2) <- paste0(rownames(adt_b2),'-adt')

rownames(rna_a1) <- paste0(rownames(rna_a1),'-rna')
rownames(rna_a3) <- paste0(rownames(rna_a3),'-rna')
rownames(rna_b3) <- paste0(rownames(rna_b3),'-rna')
rownames(rna_b2) <- paste0(rownames(rna_b2),'-rna')

assay_list = list(Age1=rbind(rna_a1,adt_a1),Age2=adt_a2,Age3=rna_a3,
                  BM1=adt_b1,BM2=rbind(rna_b2,adt_b2),BM3=rna_b3)
mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("Age1"),
               maxFeatures = 10000,
               ncomponentsReference = 20,
               ncomponentsSubset = 20)
mm <- data.frame(batch = c(rep("Age1",6316),rep("Age2",4071),rep("Age3",10639),
                           rep("BM1",9751),rep("BM2",6963),rep("BM3",11317)),
                 celltype = c(as.character(meta.list[[1]][,4]),
                              as.character(meta.list[[2]][,4]),
                              as.character(meta.list[[3]][,4]),
                              as.character(meta.list[[4]][,4]),
                              as.character(meta.list[[5]][,4]),
                              as.character(meta.list[[6]][,4])),
                 celltype1 = c(as.character(meta.list[[1]][,5]),
                               as.character(meta.list[[2]][,5]),
                               as.character(meta.list[[3]][,5]),
                               as.character(meta.list[[4]][,5]),
                               as.character(meta.list[[5]][,5]),
                               as.character(meta.list[[6]][,5])),
                 celltype2 = c(as.character(meta.list[[1]][,3]),
                               as.character(meta.list[[2]][,3]),
                               as.character(meta.list[[3]][,3]),
                               as.character(meta.list[[4]][,3]),
                               as.character(meta.list[[5]][,3]),
                               as.character(meta.list[[6]][,3])),
                 modality = c(as.character(meta.list[[1]][,6]),
                              as.character(meta.list[[2]][,6]),
                              as.character(meta.list[[3]][,6]),
                              as.character(meta.list[[4]][,6]),
                              as.character(meta.list[[5]][,6]),
                              as.character(meta.list[[6]][,6]))
)

out_dir <- './Benchmarking/Unsupervised/Ab-seq/output'
write.csv(as.matrix(stab),
          paste0(out_dir,'/stab_s4.csv'),
          row.names = FALSE)
