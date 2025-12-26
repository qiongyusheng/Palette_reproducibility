library(tools)
library(glue)
library(data.table)
library(Seurat)
library(gridExtra)
library(SingleCellExperiment)
library(parallel)
library(scran)
source('./mrtrics/imputation_helper.R')

s3_rna_df <- evalue_rna(path_real = '.',
                       path_impute = './imputation/Retina',# 预测数据的路径
                       meta_path = './meta.rds',# meta数据路径
                       group = 'cell_type__custom',# 使用的meta列名
                       name_real = 'rna.rds',# 真实数据的名字
                       name_impute = c('/Palette/LGS3OS_rna.csv',
                                       '/MIDAS/LGS3OS_rna.csv',
                                       '/Multigrate/LGS3OS_rna.csv',
                                      '/scVAEIT/LGS3OS_rna.csv'),# 各个方法的分支路径
                       name_method = c('Palette','MIDAS','Multigrate','scVAEIT'))
s3_rna_df$task <- c('Palette','MIDAS','Multigrate','scVAEIT')

d1_rna_df <- evalue_rna(path_real = '.',
                       path_impute = './imputation/Retina',# 预测数据的路径
                       meta_path = './meta.rds',# meta数据路径
                       group = 'cell_type__custom',# 使用的meta列名
                       name_real = 'rna.rds',# 真实数据的名字
                       name_impute = c('/Palette/LVG1OD_rna.csv',
                                       '/MIDAS/LVG1OD_rna.csv',
                                       '/Multigrate/LVG1OD_rna.csv',
                                      '/scVAEIT/LVG1OD_rna.csv'),# 各个方法的分支路径
                       name_method = c('Palette','MIDAS','Multigrate','scVAEIT'))
d1_rna_df$task <- c('Palette','MIDAS','Multigrate','scVAEIT')

s1_rna_df <- evalue_rna(path_real = '.',
                       path_impute = './imputation/Retina',# 预测数据的路径
                       meta_path = './meta.rds',# meta数据路径
                       group = 'cell_type__custom',# 使用的meta列名
                       name_real = 'rna.rds',# 真实数据的名字
                       name_impute = c('/Palette/LVG1OS_rna.csv',
                                       '/MIDAS/LVG1OS_rna.csv',
                                       '/Multigrate/LVG1OS_rna.csv',
                                      '/scVAEIT/LVG1OS_rna.csv'),# 各个方法的分支路径
                       name_method = c('Palette','MIDAS','Multigrate','scVAEIT'))
s1_rna_df$task <- c('Palette','MIDAS','Multigrate','scVAEIT')

rna_df <- rbind(s3_rna_df,d1_rna_df,s1_rna_df)
rna_df

d2_atac_df_1 <- evalue_atac(path_real = '.',
                       path_impute = './imputation/Retina',# 预测数据的路径
                       meta_path = './meta.rds',# meta数据路径
                       group = 'cell_type__custom',# 使用的meta列名
                       name_real = 'atac.rds',# 真实数据的名字
                       name_impute = c('Palette/LGS2OD_atac.mtx'),# 各个方法的分支路径
                       name_method = c('Palette'))

d2_atac_df_2 <- evalue_atac(path_real = '.',
                       path_impute = './imputation/Retina',# 预测数据的路径
                       meta_path = './meta.rds',# meta数据路径
                       group = 'cell_type__custom',# 使用的meta列名
                       name_real = 'atac.rds',# 真实数据的名字
                       name_impute = c('MIDAS/LGS2OD_atac.csv',
                                       'Multigrate/LGS2OD_atac.mtx'),# 各个方法的分支路径
                       name_method = c('MIDAS','Multigrate'))


d2_atac_df_3 <- evalue_atac(path_real = '.',
                       path_impute = './imputation/Retina',# 预测数据的路径
                       meta_path = './meta.rds',# meta数据路径
                       group = 'cell_type__custom',# 使用的meta列名
                       name_real = 'atac.rds',# 真实数据的名字
                       name_impute = c('scVAEIT/LGS2OD_atac.csv'),# 各个方法的分支路径
                       name_method = c('scVAEIT'),
                       feat_reduce = './imputation/Retina/scVAEIT/peak.csv')
d2_atac_df_1$task <- c('Palette')
d2_atac_df_2$task <- c('MIDAS','Multigrate')
d2_atac_df_3$task <- c('scVAEIT')
d2_atac_df <- rbind(d2_atac_df_1,d2_atac_df_2,d2_atac_df_3)

s2_atac_df_1 <- evalue_atac(path_real = './',
                       path_impute = './imputation/Retina',# 预测数据的路径
                       meta_path = '/meta.rds',# meta数据路径
                       group = 'cell_type__custom',# 使用的meta列名
                       name_real = 'atac.rds',# 真实数据的名字
                       name_impute = c('Palette/LGS2OS_atac.mtx'),# 各个方法的分支路径
                       name_method = c('Palette'))

s2_atac_df_2 <- evalue_atac(path_real = './',
                       path_impute = './imputation/Retina',# 预测数据的路径
                       meta_path = '/meta.rds',# meta数据路径
                       group = 'cell_type__custom',# 使用的meta列名
                       name_real = 'atac.rds',# 真实数据的名字
                       name_impute = c('MIDAS/LGS2OS_atac.csv',
                                       'Multigrate/LGS2OS_atac.mtx'),# 各个方法的分支路径
                       name_method = c('MIDAS','Multigrate'))


s2_atac_df_3 <- evalue_atac(path_real = './',
                       path_impute = './imputation/Retina',# 预测数据的路径
                       meta_path = '///meta.rds',# meta数据路径
                       group = 'cell_type__custom',# 使用的meta列名
                       name_real = 'atac.rds',# 真实数据的名字
                       name_impute = c('scVAEIT/LGS2OS_atac.csv'),# 各个方法的分支路径
                       name_method = c('scVAEIT'),
                       feat_reduce = './imputation/Retina/scVAEIT/peak.csv')
s2_atac_df_1$task <- c('Palette')
s2_atac_df_2$task <- c('MIDAS','Multigrate')
s2_atac_df_3$task <- c('scVAEIT')
s2_atac_df <- rbind(s2_atac_df_1,s2_atac_df_2,s2_atac_df_3)

d3_atac_df_1 <- evalue_atac(path_real = './',
                       path_impute = './imputation/Retina',# 预测数据的路径
                       meta_path = './meta.rds',# meta数据路径
                       group = 'cell_type__custom',# 使用的meta列名
                       name_real = 'atac.rds',# 真实数据的名字
                       name_impute = c('Palette/LGS3OD_atac.mtx'),# 各个方法的分支路径
                       name_method = c('Palette'))

d3_atac_df_2 <- evalue_atac(path_real = './',
                       path_impute = './imputation/Retina',# 预测数据的路径
                       meta_path = './meta.rds',# meta数据路径
                       group = 'cell_type__custom',# 使用的meta列名
                       name_real = 'atac.rds',# 真实数据的名字
                       name_impute = c('MIDAS/LGS3OD_atac.csv',
                                       'Multigrate/LGS3OD_atac.mtx'),# 各个方法的分支路径
                       name_method = c('MIDAS','Multigrate'))


d3_atac_df_3 <- evalue_atac(path_real = './',
                       path_impute = './imputation/Retina',# 预测数据的路径
                       meta_path = './meta.rds',# meta数据路径
                       group = 'cell_type__custom',# 使用的meta列名
                       name_real = 'atac.rds',# 真实数据的名字
                       name_impute = c('scVAEIT/LGS3OD_atac.csv'),# 各个方法的分支路径
                       name_method = c('scVAEIT'),
                       feat_reduce = './imputation/Retina/scVAEIT/peak.csv')
d3_atac_df_1$task <- c('Palette')
d3_atac_df_2$task <- c('MIDAS','Multigrate')
d3_atac_df_3$task <- c('scVAEIT')
d3_atac_df <- rbind(d3_atac_df_1,d3_atac_df_2,d3_atac_df_3)

atac_df <- rbind(d2_atac_df,s2_atac_df,d3_atac_df)
atac_df

res <- list(RNA = rna_df,
           ATAC = atac_df)
saveRDS(res,'./imputation/Retina/metrics.rds')
