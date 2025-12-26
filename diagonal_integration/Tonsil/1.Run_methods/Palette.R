source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

setwd('./')

rna <- readRDS('./RNA/rna.rds')
codex <- readRDS('./CODEX/codex.rds')
rna_co <- readRDS('./RNA/co_clean_smooth_q.rds')
codex_co <- readRDS('./CODEX/co_clean_smooth_q.rds')

colnames(codex_co) <- colnames(codex)

musa_obj <- Create.Musa.Object(data.list = list(rna,
                                                rna_co,
                                                codex,
                                                codex_co), 
                               samples = c('B1','B1','B2','B2'), 
                               modals = c("rna","co","adt","co"),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna",'co','adt'), 
                               sparce = TRUE)

musa_obj@Assays[[1]]@data <- musa_obj@Assays[[1]]@Raw
musa_obj@Assays[[2]]@data <- musa_obj@Assays[[2]]@Raw
musa_obj@Assays[[3]]@data <- musa_obj@Assays[[3]]@Raw
musa_obj@Assays[[4]]@data <- musa_obj@Assays[[4]]@Raw

musa_obj <- Add.HVFs(musa_obj,modal = "rna")
musa_obj <- Add.HVFs(musa_obj, modal = "adt")
musa_obj <- Add.HVFs(musa_obj, modal = "co")

musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("co"),
                         sample = NULL,
                         dim.reduce = TRUE,# 是否在降维矩阵上计算
                         dims = list(1:30),# 使用哪些维数
                         method = c("PCA"),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)

musa_obj <- Find.Subspace3(musa_obj,
                          modal = c("co"),
                          lambda = list(0.9),
                          alpha = list(0),
                         sub.dims = list(1:20),
                         cos.dims = list(1:20),
                          Angle.var = c(30),
                         max.Angle = c(60))

musa_obj <- IsolatedModalEmbedding(musa_obj,modal = c("adt",'rna'),
                                   dims = list(1:30,1:50))

rna <- readMM('/home/server/sqy/data/mosaic/RAW/Human_Tonsil_RNA_CODEX_MaxFuse/RNA/tonsil_rna_0510.txt')
gene_name <- fread('/home/server/sqy/data/mosaic/RAW/Human_Tonsil_RNA_CODEX_MaxFuse/RNA/tonsil_rna_0510_names.csv',
                   data.table = F,header = T)
rna_meta <- fread('/home/server/sqy/data/mosaic/RAW/Human_Tonsil_RNA_CODEX_MaxFuse/RNA/tonsil_rna_0510_meta.csv',
                  data.table = F)
rna <- as(t(rna),'dgCMatrix')
rownames(rna) <- gene_name[,2]
colnames(rna) <- rna_meta[,1]
musa_obj@Assays[[1]]@Raw <- rna
rna <- Seurat::NormalizeData(rna)
musa_obj@Assays[[1]]@data <- rna
musa_obj <- Add.HVFs(musa_obj,modal = "rna")

a <- Run.Palette2(musa_obj,nn = 2,lambda = 0.5,nDims = 20, modal.norm = F, origin.infer = T)

infer_rna <- a@Assays$rna_B2@data

out_dir <- './'

write.csv(as.matrix(a@Int.result[["bind"]]),
          paste0(out_dir,'Palette.csv'),
          row.names = FALSE)

writeMM(infer_rna,
        './infer_rna.mtx')
write.csv(rownames(infer_rna),
          './gene.csv')
write.csv(colnames(infer_rna),
          './barcode.csv')