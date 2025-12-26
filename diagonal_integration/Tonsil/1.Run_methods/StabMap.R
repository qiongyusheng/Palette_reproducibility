library(bindSC)
library(Seurat)
library(data.table)
library(Matrix)
library(matrixStats)
library(StabMap)
library(scran)
library(patchwork)

setwd('./')
rna <- as(as.matrix(readRDS('./RNA/rna.rds')),'dgCMatrix')
# cytof <- as(as.matrix(readRDS('./cytof/cytof.rds')),'dgCMatrix')
rna_co <- as(as.matrix(readRDS('./RNA/co_all.rds')),'dgCMatrix')
# codex_co <- as.matrix(readRDS('./CODEX/co_all.rds'))
codex <- as.matrix(readRDS('./CODEX/codex.rds'))
co_feat <- readRDS('./co_feat.rds')

remain_gene <- setdiff(rownames(rna),co_feat[,1])
rna <- rna[remain_gene,]
# rna <- rbind(rna,rna_co)
rownames(rna_co) <- co_feat[,1]
rna <- rbind(rna,rna_co)

codex_co <- codex[co_feat[,2],]
remain_pro <- setdiff(rownames(codex),co_feat[,2])
codex <- codex[remain_pro,]
rownames(codex_co) <- co_feat[,1]
codex <- rbind(codex,codex_co)

assay_list = list(RNA=rna,CODEX=codex)
mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("CODEX"),
               maxFeatures = 10000,
              ncomponentsReference = 20,
              ncomponentsSubset = 20)
out <- rbind(stab[178920:nrow(stab),],stab[1:178919,])
out_dir <- './'
write.csv(as.matrix(out),
          paste0(out_dir,'/stab.csv'),
          row.names = FALSE)