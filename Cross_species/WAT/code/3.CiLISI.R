library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics_R.R')

data_dir <- './'
meta_dir <- './'
out_dir <- './'


method = c('Multigrate_homo','Multigrate_non',
           'scMoMaT_homo','scMoMaT_non','midas_homo','midas_non','scVAEIT_homo',
           'scVAEIT_non',"stab_non","stab_homo",'scpoli','SIGNAL','Seurat',
           'scanvi','Harmony','fastMNN',"Palette_homo","Palette_non")

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0(data_dir,method[i],'.csv'),data.table = F,header = T)
}
meta <- readRDS(paste0(meta_dir,'meta.rds'))

out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'Species','CellType',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c('Multigrate_homo',
                  'Multigrate_non',
                  'scMoMaT_homo',
                  'scMoMaT_non',
                 'MIDAS_homo',
                  'MIDAS_non',
                 'scVAEIT_homo',
                  'scVAEIT_non',
                 "StabMap_non",
                 "StabMap_homo",
                 'scPoli',
                 'SIGNAL',
                 'Seurat',
                 'scANVI',
                 'Harmony',
                  'fastMNN',"Palette_homo",
                  "Palette_non")
saveRDS(out,'./lisi_score.rds')
