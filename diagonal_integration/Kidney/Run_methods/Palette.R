source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

setwd('./')

rna <- readRDS('./RNA/pca.rds')
atac <- readRDS('./ATAC/lsi.rds')
rna_co <- readRDS('./RNA/co_n_q.rds')
atac_co <- readRDS('./ATAC/co_n_q.rds')

musa_obj <- Create.Musa.Object(data.list = list(t(rna),rna_co,t(atac),atac_co), 
                               samples = c('B1','B1','B2','B2'), 
                               modals = c("rna","co","atac","co"),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna",'co','atac'), 
                               sparce = TRUE)

options(future.globals.maxSize = 80000 * 1024^2)
musa_obj@Assays[[1]]@data <- musa_obj@Assays[[1]]@Raw
musa_obj@Assays[[2]]@data <- musa_obj@Assays[[2]]@Raw
musa_obj@Assays[[3]]@data <- musa_obj@Assays[[3]]@Raw
musa_obj@Assays[[4]]@data <- musa_obj@Assays[[4]]@Raw

musa_obj <- Add.HVFs(musa_obj, modal = "atac")
musa_obj <- Add.HVFs(musa_obj,modal = "rna")
musa_obj <- Add.HVFs(musa_obj, modal = "co")

musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("co"),
                         dims = list(1:50),
                         method = c("PCA"))

musa_obj@Assays[[3]]@embedding <- as.matrix(musa_obj@Assays[[3]]@Raw[1:30,])
musa_obj@Assays[[1]]@embedding <- as.matrix(musa_obj@Assays[[1]]@Raw[1:30,])
musa_obj <- Find.Subspace3(musa_obj,
                          modal = c("co"),
                          lambda = list(0.9),
                         sub.dims = list(1:30),
                         cos.dims = list(1:50),
                          Angle.var = c(30),
                         max.Angle = c(60))

a <- Run.Palette2(musa_obj, modal.norm = F)

meta <- readRDS('./meta.rds')

out_dir <- './'
write.csv(as.matrix(a@Int.result[["bind"]]),
          paste0(out_dir,'Palette.csv'),
          row.names = FALSE)