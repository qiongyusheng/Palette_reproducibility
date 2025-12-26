library(StabMap)
library(scran)
library(Seurat)
library(patchwork)
setwd("./")
path <- c('./PBMC_ASAP/control_asap','./PBMC_ASAP/stim_asap',
          './PBMC_CITE/control_cite','./PBMC_CITE/stim_cite')
atac.list <- list()
adt.list <- list()
rna.list <- list()
meta.list <- list()
for(i in 1:4){
    meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
    adt.list[[i]] <- readRDS(paste0(path[i],'/adt.rds'))
    if(i %in% c(1,2)){
        atac.list[[i]] <- readRDS(paste0(path[i],'/atac.rds'))
    }
    if(i %in% c(3,4)){
        rna.list[[i-2]] <- readRDS(paste0(path[i],'/rna.rds'))
    }
}

for(i in 1:4){
  adt.list[[i]] <- NormalizeData(adt.list[[i]])
}
for(i in 1:2){
  rna.list[[i]] <- NormalizeData(rna.list[[i]])
  atac.list[[i]] <- NormalizeData(atac.list[[i]])
}

# adt for asap and cite
adt1 <- adt.list[[1]]
adt2 <- adt.list[[2]]
adt3 <- adt.list[[3]]
adt4 <- adt.list[[4]]

# atac for asap
atac1 <- atac.list[[1]]
atac2 <- atac.list[[2]]

# rna for cite
rnac3 <- rna.list[[1]]
rnac4 <- rna.list[[2]]

rownames(rnac3) <- rownames(rnac4) <- paste0(rownames(rnac3),'-rna')
rownames(adt1) <- rownames(adt2) <- rownames(adt3) <- rownames(adt4) <- paste0(rownames(adt1),'-adt')
rownames(atac1) <- rownames(atac2) <- paste0(rownames(atac1),'-atac')


assay_list = list(asap_control=rbind(adt1,atac1),
                  asap_stim=rbind(adt2,atac2),
                  cite_control=rbind(adt3,rnac3),
                  cite_stim=rbind(adt4,rnac4))
mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("cite_control"),
               maxFeatures = 1000000,
              ncomponentsReference = 20,
              ncomponentsSubset = 20)

Tname <- c(colnames(adt1),colnames(adt2),colnames(adt3),colnames(adt4))
stab <- stab[Tname,]
write.csv(as.matrix(stab),
          "./Cross_condition_PBMC/output/Embedding/stab.csv",
          row.names = FALSE)
