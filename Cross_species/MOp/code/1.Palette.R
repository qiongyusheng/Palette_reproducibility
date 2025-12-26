source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
options(future.globals.maxSize = 500 * 1024^3)

setwd('./cross_species')
Run_Multi_SVD <- function(x, n.dims = 50, thresh = 10000L){
  if(nrow(x[[1]]) < thresh){
    
    d <- lapply(x, function(a){
      a <- tcrossprod(a)/(ncol(a)-1)
      a
    })
    
    d <- Reduce("+", d)
    V <- RSpectra::svds(d, k = n.dims)[["u"]] %>% as.matrix()
    r <- lapply(x, function(b){
      crossprod(V, b)
    })
    return(r)
    
  }else{
    d <- lapply(x, function(a){
      a <- a/sqrt(ncol(a)-1)
      a
    })
    d <- Reduce(cbind, d)
    V <- RSpectra::svds(d, k = n.dims)[["u"]] %>% as.matrix()
    r <- lapply(x, function(b){
      crossprod(V, b)
    })
    return(r)
  }
}

# human
# 1.SNARE (RNA+ATAC)
h_s_atac <- as(readRDS('./human/SNARE_Lein/atac.rds'),'dgCMatrix')
h_s_rna <- as(as.matrix(readRDS('./human/SNARE_Lein/co.rds')),'dgCMatrix')
h_s_m <- readRDS('./human/SNARE_Lein/m.rds')

# 2.10XV3 RNA
h_10xv3_rna <- as(as.matrix(readRDS('./human/10XV3_Lein/co.rds')),'dgCMatrix')
h_10xv3_m <- readRDS('./human/10XV3_Lein/m.rds')

# 3. 10XMultiome ATAC
h_10xmulti_atac <- as(readRDS('./human/Multiome_Zemke/d1_atac.rds'),'dgCMatrix')
h_10xmulti_m <- readRDS('./human/Multiome_Zemke/m.rds')
rownames(h_10xmulti_atac) <- rownames(h_s_atac)

# mouse
# 1.10XMultiome (RNA+ATAC)
m_10xmulti_rna <- as(as.matrix(readRDS('./mouse/10XMultiome_Zemke/co.rds')),'dgCMatrix')
m_10xmulti_atac <- as(as.matrix(readRDS('./mouse/10XMultiome_Zemke/atac.rds')),'dgCMatrix')
m_10xmulti_m <- readRDS('./mouse/10XMultiome_Zemke/m.rds')

# 2. 10XV3 RNA
m_10xv3_rna <- as(as.matrix(readRDS('./mouse/10XV3_Lein/co.rds')),'dgCMatrix')
m_10xv3_m <- readRDS('./mouse/10XV3_Lein/m.rds')

# 3. 10X ATAC
m_10x_atac <- as(as.matrix(readRDS('./mouse/10X/atac.rds')),'dgCMatrix')
m_10x_m <- readRDS('./mouse/10X/m.rds')

# Marmoset
# 1.SNARE (RNA+ATAC)
mar_snare_rna <- as(as.matrix(readRDS('./Marmoset/SNARE_Lein/co.rds')),'dgCMatrix')
mar_snare_m <- readRDS('./Marmoset/SNARE_Lein/m.rds')

# 2.10XV3 RNA
mar_10xv3_rna <- as(as.matrix(readRDS('./Marmoset/10XV3_Lein/co.rds')),'dgCMatrix')
mar_10xv3_m <- readRDS('./Marmoset/10XV3_Lein/m.rds')

gc()

musa_obj <- Create.Musa.Object(data.list = list(h_s_atac,h_s_rna,
                                                h_10xv3_rna,
                                                h_10xmulti_atac,
                                                m_10xmulti_rna,m_10xmulti_atac,
                                                m_10xv3_rna,
                                                m_10x_atac,
                                                mar_snare_rna,
                                                mar_10xv3_rna), 
                               samples = c('human_snare','human_snare',
                                           'human_10xv3',
                                           'human_10xmulti',
                                           'mouse_10xmulti','mouse_10xmulti',
                                           'mouse_10xv3',
                                           'mouse_10xatac',
                                           'mar_snare',
                                           'mar_10xv3'), 
                               modals = c('h-atac','rna',
                                          'rna',
                                          'h-atac',
                                          'rna','m-atac',
                                          'rna',
                                          'm-atac',
                                          'rna',
                                          'rna'),
                               filter = F,
                               min.cells = c(0,0,0,0),  
                               min.features = c(1,1,1,1),
                               modal.order = c("rna","h-atac",'m-atac'), 
                               sparce = TRUE)
rm(h_s_atac,h_s_rna,
   h_10xv3_rna,
   h_10xmulti_atac,
   m_10xmulti_rna,m_10xmulti_atac,
   m_10xv3_rna,
   m_10x_atac,
   mar_snare_rna,
   mar_10xv3_rna)
gc()

options(future.globals.maxSize = 50 * 1024^3)
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "rna",
                           normal.method = "LogNormalize")
musa_obj <- IDFLog_try(musa_obj,modal = "h-atac")
musa_obj <- IDFLog_try(musa_obj,modal = "m-atac")

musa_obj <- Add.HVFs(musa_obj,modal = "rna")
musa_obj <- Add.HVFs(musa_obj, modal = "h-atac")
musa_obj <- Add.HVFs(musa_obj, modal = "m-atac")

musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("rna",'h-atac','m-atac'),
                         sample = NULL,
                         dim.reduce = TRUE,# 是否在降维矩阵上计算
                         dims = list(1:50,1:50,1:50),# 使用哪些维数
                         method = c("PCA",'LSI','LSI'),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'h-atac','m-atac'),
                           lambda = list(0.8,0.5,0.5),
                           joint = F,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30,1:30),
                           cos.dims = list(1:50,1:50,1:50),
                           Angle.var = c(15,20,20),
                           max.Angle = c(50,50,50),
                           pre.reduce = TRUE,
                           pre_dims = list(1:100))

a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,emb_L2 =T)

emb <- a@Int.result$bind
meta <- rbind(h_s_m,h_10xv3_m,h_10xmulti_m,
              m_10xmulti_m,m_10xv3_m,m_10x_m,
              mar_snare_m,mar_10xv3_m)
emb <- emb[rownames(meta),]
out_dir <- './'
write.csv(as.matrix(emb),
          paste0(out_dir,'/Palette.csv'),
          row.names = FALSE)
