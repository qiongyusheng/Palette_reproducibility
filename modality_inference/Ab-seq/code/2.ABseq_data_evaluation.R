library(tools)
library(glue)
library(data.table)
library(Seurat)
library(gridExtra)
library(SingleCellExperiment)
library(parallel)
library(scran)
source('./metrics/imputation_helper.R')

a2_rna_df <- evalue_rna(path_real = './Young_Age/Age2/',
                       path_impute = './ABseq/predict',# 预测数据的路径
                       meta_path = './Young_Age/Age2/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'rna.rds',# 真实数据的名字
                       name_impute = c('/palette/Age2rna_k10.csv',
                                       '/midas/Age2rna.csv',
                                       '/multigrate/Age2rna.csv',
                                      '/scvaeit/Age2rna.csv'),# 各个方法的分支路径
                       name_method = c('Palette','MIDAS','Multigrate','scVAEIT'))
a2_rna_df$task <- c('Palette','MIDAS','Multigrate','scVAEIT')

m3_rna_df <- evalue_rna(path_real = './Young_Age/BM3/',
                       path_impute = './ABseq/predict',# 预测数据的路径
                       meta_path = './Young_Age/BM3/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'rna.rds',# 真实数据的名字
                       name_impute = c('/palette/BM3rna_k10.csv',
                                       '/midas/BM3rna.csv',
                                       '/multigrate/BM3rna.csv',
                                      '/scvaeit/BM3rna.csv'),# 各个方法的分支路径
                       name_method = c('Palette','MIDAS','Multigrate','scVAEIT'))
m3_rna_df$task <- c('Palette','MIDAS','Multigrate','scVAEIT')

rna_df <- rbind(a2_rna_df,m3_rna_df)

a3_adt_df <- evalue_adt(path_real = './Young_Age/Age3/',
                       path_impute = './ABseq/predict',# 预测数据的路径
                       meta_path = './Young_Age/Age3/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'adt.rds',# 真实数据的名字
                       name_impute = c('/palette/Age3adt_k10.csv',
                                       '/midas/Age3adt.csv',
                                       '/multigrate/Age3adt.csv',
                                      '/scvaeit/Age3adt.csv'),# 各个方法的分支路径
                       name_method = c('Palette','MIDAS','Multigrate','scVAEIT'))

m2_adt_df <- evalue_adt(path_real = './Young_Age/BM2/',
                       path_impute = './ABseq/predict',# 预测数据的路径
                       meta_path = './Young_Age/BM2/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'adt.rds',# 真实数据的名字
                       name_impute = c('/palette/BM2adt_k10.csv',
                                       '/midas/BM2adt.csv',
                                       '/multigrate/BM2adt.csv',
                                      '/scvaeit/BM2adt.csv'),# 各个方法的分支路径
                       name_method = c('Palette','MIDAS','Multigrate','scVAEIT'))

a3_adt_df$task <- c('Palette','MIDAS','Multigrate','scVAEIT')
m2_adt_df$task <- c('Palette','MIDAS','Multigrate','scVAEIT')

adt_df <- rbind(a3_adt_df,m2_adt_df)

res <- list(RNA = rna_df,
           ADT = adt_df)
saveRDS(res,'./ABseq/metrics.rds')
