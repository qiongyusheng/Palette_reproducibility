source("./midas/preprocess/utils.R")
# setwd("/root/workspace/code/midas/")
library(RColorBrewer)

impute_out <- './'

data_name <- c('B3_adt','B3_rna','B4_adt','B4_atac','B5_atac','B5_rna')

path <- c('/midas/result/TEA_s1_impu/e0/default/predict/sp_latest/subset_1/x_impt/adt/',
         '/midas/result/TEA_s1_impu/e0/default/predict/sp_latest/subset_1/x_impt/rna/',
         '/midas/result/TEA_s1_impu/e0/default/predict/sp_latest/subset_2/x_impt/adt/',
         '/midas/result/TEA_s1_impu/e0/default/predict/sp_latest/subset_2/x_impt/atac/',
         '/midas/result/TEA_s1_impu/e0/default/predict/sp_latest/subset_3/x_impt/atac/',
         '/midas/result/TEA_s1_impu/e0/default/predict/sp_latest/subset_3/x_impt/rna/')



for(i in 1:6){
    miss_dir <- path[i]
    fnames1 <- dir(path = miss_dir, pattern = ".csv$")
    fnames1 <- str_sort(fnames1, decreasing = F)
    z_subset_list <- list()
    N <- length(fnames1)
    
    for (n in seq_along(fnames1)) {
        message(paste0("Loading impute data ", i, "/8", ", File ", n, "/", N))
        z_subset_list[[n]] <- read.csv(file.path(miss_dir, fnames1[n]), header = F)
      }  
    
    tmp <- bind_rows(z_subset_list)
    print(dim(tmp))
    write.csv(tmp,paste0(impute_out,data_name[i],'.csv'))
              }

