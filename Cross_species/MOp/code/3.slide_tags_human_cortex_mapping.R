library(Seurat)
library(dplyr)
library(data.table)

data <- fread("./slide_tags_human_cortex/humancortex_expression.csv.gz",
              data.table = F)
rownames(data) <- data[,1]
data <- data[,-1]
data <- as(as.matrix(data),'dgCMatrix')

meta <- fread('./slide_tags_human_cortex/humancortex_metadata.csv',
              data.table = F)
meta <- meta[-1,]
rownames(meta) <- meta[,1]

coord <- fread('./slide_tags_human_cortex/humancortex_spatial.csv',
               data.table = F)
coord <- coord[-1,]
rownames(coord) <- coord[,1]
all.equal(rownames(meta),rownames(coord))
meta$X <- coord$X
meta$Y <- coord$Y

setwd('./cross_species')
# SNARE 
h_snare <- readRDS('./human/SNARE_Lein/rna.rds')
h_snare_m <- readRDS('./human/SNARE_Lein/meta.rds')
# 10Xv3
h_10x <- readRDS('./human/10XV3_Lein/rna1.rds')
h_10x_m <- readRDS('./human/10XV3_Lein/meta1.rds')
h_snare <- as(as.matrix(h_snare),'dgCMatrix')
h_10x <- as(as.matrix(h_10x),'dgCMatrix')
h_co_gene <- intersect(rownames(h_snare),rownames(h_10x))
h_snare <- h_snare[h_co_gene,]
h_10x <- h_10x[h_co_gene,]

h_co_gene <- intersect(rownames(h_snare),rownames(data))
length(h_co_gene)

h_snare <- h_snare[h_co_gene,]
h_10x <- h_10x[h_co_gene,]
data <- data[h_co_gene,]

h_snare_obj <- CreateSeuratObject(h_snare)
h_snare_obj <- NormalizeData(h_snare_obj)
h_snare_obj <- FindVariableFeatures(h_snare_obj,nfeatures = 3000)

h_10x_obj <- CreateSeuratObject(h_10x)
h_10x_obj <- NormalizeData(h_10x_obj)
h_10x_obj <- FindVariableFeatures(h_10x_obj,nfeatures = 3000)

obj <- CreateSeuratObject(data)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj,nfeatures = 3000)

h_hvg <- Seurat::SelectIntegrationFeatures(list(h_snare_obj,h_10x_obj,obj), nfeatures = 3000)
rm(h_snare_obj,h_10x_obj,obj)

h_snare <- h_snare[h_hvg,]
h_10x <- h_10x[h_hvg,]
data <- data[h_hvg,]

source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

latent_data <- t(as.matrix(fread('./Cross_species/MOp/output/Palette.csv',data.table = F)))
meta <- readRDS('/home/server/sqy/data/mosaic/result/cross_species/meta.rds')
colnames(latent_data) <- rownames(meta)

options(future.globals.maxSize = 50 * 1024^3)
res <- Palette_map(list(data,h_snare,h_10x),
                   modal = c('rna','rna','rna'),
                  batch = c('slide_tags','snare','10X'),
                  q_or_batch = c('query',rep('ref',2)),
                  modal_order = c('rna'),
                  method = 'fastMNN',
                   method_out.npcs = 40,
                  normalize_method = c('LogNormalize'),
                  latent = latent_data,Bi_sPCA = T)

result_fastMNN <- run_fastMNN(list(res[[1]],latent_data),
                      cell_name = c(colnames(res[[1]]),colnames(latent_data)),
                      is.normalize = F,out.npcs = NA)
meta$QR <- 'reference'

meta_q <- fread('./slide_tags_human_cortex/humancortex_metadata.csv',
              data.table = F)
meta_q <- meta_q[-1,]
rownames(meta_q) <- meta_q[,1]
meta_q$cellID <- rownames(meta_q)
meta_q$anno <- meta_q$cluster
meta_q$species <- 'Human'
meta_q$Modal <- 'RNA'
meta_q$Batch <- 'slide-tags'
meta_q$QR <- 'query'

meta_q <- meta_q[,c('cellID', 'anno', 'species', 'Modal', 'Batch', 'QR')]
mm <- rbind(meta[,c('cellID', 'anno', 'species', 'Modal', 'Batch', 'QR')],meta_q)

result_fastMNN <- result_fastMNN[,rownames(mm)]

Trans.Knowledge <- function(ref,
                            query,
                            ref.label = NULL,
                            k = 5){
  
  pred_label <- class::knn(ref,
                           query,
                           as.character(ref.label[,1]),
                           k = k, prob = TRUE) %>% as.character()
  
  out <- data.frame(predict_label = pred_label)
  
  return(out)
  
}

pred_label <- Trans.Knowledge(t(result_fastMNN[,rownames(meta)]),
                                         t(result_fastMNN[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta),'anno']))

mm$pred_label <- c(meta$anno,pred_label[,1])

motif1 <- readRDS('./cross_species/human/SNARE_Lein/motif.rds')
motif2 <- readRDS('./cross_species/human/Multiome_Zemke/d1_motif.rds')

motif <- cbind(motif1,motif2)
rm(motif1,motif2)

co_names <- intersect(colnames(motif),colnames(result_fastMNN))
ref_sub <- result_fastMNN[,co_names]

motif <- motif[,co_names]

query <- result_fastMNN[,rownames(meta_q)]

library(BiocNeighbors)

anchors <- queryKNN(t(ref_sub),t(query),k=5)

dist <- anchors[[2]]
dist <- exp(-dist)
dd <- rowSums(dist)

dist <- diag(1/dd) %*% dist

pred_motif <- matrix(0,nrow = nrow(motif),ncol = ncol(query))
rownames(pred_motif) <- rownames(motif)
colnames(pred_motif) <- colnames(query)

for(i in 1:nrow(anchors[[1]])){
    tmp <- motif[,anchors[[1]][i,]] %*% t(dist[i,,drop = F])
    pred_motif[,i] <- tmp
}

out_dir <- './Cross_species/MOp/output/ref_query_int'

write.csv(as.matrix(t(result_fastMNN)),
          paste0(out_dir,'/result_fastMNN.csv'),
          row.names = FALSE)
saveRDS(mm,'./Cross_species/MOp/output/ref_query_int/meta.rds')
saveRDS(pred_motif,'./pred_motif.rds')
