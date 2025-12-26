library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Supervised_using_modified_label/Label_shuffled/output/Retina')
data_name <- c('Palette_s1','Palette_s2','Palette','Palette_s3','Palette_s4')
meta_name <- c('meta_s1','meta_s2','meta','meta_s3','meta_s4')
mod <- c('_shuffled5','_shuffled10','_shuffled20')
lisi_vec <- c()
out_name <- c()
for(i in 1:length(data_name)){
  meta <- readRDS(paste0('./',meta_name[i],'.rds'))
  for(j in 1:length(mod)){
    res <- fread(paste0('./',data_name[i],mod[j],'.csv'),data.table = F,header = T)
    lisi_vec <- c(lisi_vec,
                  lisi <- CiLISI_sample(res,meta,'modality','celltype',n_repeat = 4))
    out_name <- c(out_name,paste0(data_name[i],mod[j]))
  }
}
out <- data.frame(CiLISI = rep(0,length(data_name)*length(mod)))
rownames(out) <- out_name
out[,'CiLISI'] <- lisi_vec

out_dir <- './Benchmarking/Supervised_using_modified_label/Label_shuffled/output/Retina'
saveRDS(out,paste0(out_dir,'/lisi_score.rds'))
