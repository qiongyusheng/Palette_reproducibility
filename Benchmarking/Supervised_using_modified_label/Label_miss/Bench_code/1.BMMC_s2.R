library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Supervised_using_modified_label/Label_miss/output/BMMC_s2')
data_name <- c('Palette_34','Palette_35','Palette_36','Palette_24','Palette_26')
meta_name <- c('meta_34','meta_35','meta_36','meta_24','meta_26')
mod <- c('_miss5','_miss10','_miss20')
lisi_vec <- c()
out_name <- c()
for(i in 1:length(data_name)){
  meta <- readRDS(paste0('./',meta_name[i],'.rds'))
  for(j in 1:length(mod)){
    res <- fread(paste0('./',data_name[i],mod[j],'.csv'),data.table = F,header = T)
    lisi_vec <- c(lisi_vec,
                  lisi <- CiLISI_sample(res,meta,'modality','celltype.l2',n_repeat = 4))
    out_name <- c(out_name,paste0(data_name[i],mod[j]))
  }
}
out <- data.frame(CiLISI = rep(0,length(data_name)*length(mod)))
rownames(out) <- out_name
out[,'CiLISI'] <- lisi_vec

out_dir <- './Benchmarking/Supervised_using_modified_label/Label_miss/output/BMMC_s2'
saveRDS(out,paste0(out_dir,'/lisi_score.rds'))
