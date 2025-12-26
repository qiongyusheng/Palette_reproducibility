source("./midas/preprocess/utils.R")
library(RColorBrewer)
library(Matrix)

impute_out <- './'

data_name <- c('LGS2OD_atac','LGS2OS_atac','LGS3OD_atac','LGS3OS_rna','LVG1OD_rna','LVG1OS_rna')

path <- c('./midas/result/retina_impu/e0/default/predict/sp_latest/subset_2/x_impt/atac/',
         './midas/result/retina_impu/e0/default/predict/sp_latest/subset_3/x_impt/atac/',
         './midas/result/retina_impu/e0/default/predict/sp_latest/subset_4/x_impt/atac/',
         './midas/result/retina_impu/e0/default/predict/sp_latest/subset_5/x_impt/rna/',
         './midas/result/retina_impu/e0/default/predict/sp_latest/subset_6/x_impt/rna/',
         './midas/result/retina_impu/e0/default/predict/sp_latest/subset_7/x_impt/rna/')


for(i in 4:6){
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

library(data.table)
for(i in 1:3){
    miss_dir <- path[i]
    fnames1 <- dir(path = miss_dir, pattern = ".csv$")
    fnames1 <- str_sort(fnames1, decreasing = F)
    z_subset_list <- list()
    N <- length(fnames1)
    for (n in seq_along(fnames1)) {
        message(paste0("Loading impute data ", i, "/8", ", File ", n, "/", N))
        all_mat <- fread(file.path(miss_dir, fnames1[n]),data.table = F,header = F)
        all_mat <- all_mat
        z_subset_list[[n]] <- all_mat
      } 
    
    tmp <- bind_rows(z_subset_list)
    print(dim(tmp))
    write.csv(tmp,paste0(impute_out,data_name[i],'.csv'))
}
