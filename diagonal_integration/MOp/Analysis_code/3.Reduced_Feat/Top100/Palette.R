source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

setwd('./MERFISH_ATAC')

lsi <- readRDS('./ATAC/SVD30.rds')
co_atac <- readRDS('./ATAC/co_smooth_qn100_RF.rds')
meta_atac <- readRDS('./ATAC/meta.rds')

merf <- readRDS('./MERFISH/pca30.rds')
co_merf <- readRDS('./MERFISH/co_qn100_RF.rds')
meta_merf <- readRDS('./MERFISH/meta.rds')
options(future.globals.maxSize = 30 * 1024^3)
musa_obj <- Create.Musa.Object(data.list = list(t(merf),
                                                co_merf,co_atac,
                                                t(lsi)), 
                               samples = c('b1','b1','b2','b2'), 
                               modals = c('rna','co','co','atac'),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna",'co','atac'), 
                               sparce = TRUE)

musa_obj@Assays[[1]]@data <- musa_obj@Assays[[1]]@Raw
musa_obj@Assays[[2]]@data <- musa_obj@Assays[[2]]@Raw
musa_obj@Assays[[3]]@data <- musa_obj@Assays[[3]]@Raw
musa_obj@Assays[[4]]@data <- musa_obj@Assays[[4]]@Raw

musa_obj <- Add.HVFs(musa_obj,modal = "rna")
musa_obj <- Add.HVFs(musa_obj,modal = "co")
musa_obj <- Add.HVFs(musa_obj,modal = "atac")

musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("co"),
                         sample = NULL,
                         dim.reduce = TRUE,
                         dims = list(1:50),
                         method = c("PCA"),
                         joint = F,
                         resolution = 1,
                         nn.k = 20,
                         nn.method = "annoy",
                         annoy.metric = "euclidean",
                         verbose = TRUE)

musa_obj <- Find.Subspace3(musa_obj,
                          modal = c("co"),
                          lambda = list(0.9),
                          alpha = list(0),
                         sub.dims = list(1:15),
                         cos.dims = list(1:50),
                          Angle.var = c(30),
                         max.Angle = c(60))

musa_obj@Assays[[4]]@embedding <- as.matrix(musa_obj@Assays[[4]]@Raw)
musa_obj@Assays[[1]]@embedding <- as.matrix(musa_obj@Assays[[1]]@Raw)

a <- Run.Palette2(musa_obj,nn = 2,lambda = 0.5,nDims = 20,modal.norm = F)
# save data
out_dir <- './Reduced_feat/RF100/'
write.csv(as.matrix(a@Int.result[["bind"]]),
          paste0(out_dir,'Palette.csv'),
          row.names = FALSE)


