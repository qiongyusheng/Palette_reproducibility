library(tools)
library(glue)
library(data.table)
library(Seurat)
library(gridExtra)
library(SingleCellExperiment)
library(parallel)
library(scran)
source('./metrics/imputation_helper.R')

B3_rna_df_1 <- evalue_rna(path_real = './PBMC_TEA/B3/',
                       path_impute = './imputation/TEA_s1',# 预测数据的路径
                       meta_path = './PBMC_TEA/B3/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'rna.rds',# 真实数据的名字
                       name_impute = c('/Palette/B3_rna.csv',
                                       '/MIDAS/B3_rna.csv',
                                       '/Multigrate/B3_rna.csv'),# 各个方法的分支路径
                       name_method = c('Palette','MIDAS','Multigrate'))
B3_rna_df_1$task <- c('Palette','MIDAS','Multigrate')

B3_rna_df_2 <- evalue_rna(path_real = './PBMC_TEA/B3/',
                       path_impute = './imputation/TEA_s1',# 预测数据的路径
                       meta_path = './PBMC_TEA/B3/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'rna.rds',# 真实数据的名字
                       name_impute = c('/scVAEIT/B3RNA.csv'),# 各个方法的分支路径
                       name_method = c('scVAEIT'),
                      feat_reduce = './imputation/TEA_s1/scVAEIT/gene.csv')
B3_rna_df_2$task <- c('scVAEIT')

B5_rna_df_1 <- evalue_rna(path_real = './PBMC_TEA/B5/',
                       path_impute = './imputation/TEA_s1',# 预测数据的路径
                       meta_path = './PBMC_TEA/B5/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'rna.rds',# 真实数据的名字
                       name_impute = c('/Palette/B5_rna.csv',
                                       '/MIDAS/B5_rna.csv',
                                       '/Multigrate/B5_rna.csv'),# 各个方法的分支路径
                       name_method = c('Palette','MIDAS','Multigrate'))
B5_rna_df_1$task <- c('Palette','MIDAS','Multigrate')

B5_rna_df_2 <- evalue_rna(path_real = './PBMC_TEA/B5/',
                       path_impute = './imputation/TEA_s1',# 预测数据的路径
                       meta_path = './PBMC_TEA/B5/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'rna.rds',# 真实数据的名字
                       name_impute = c('/scVAEIT/B5RNA.csv'),# 各个方法的分支路径
                       name_method = c('scVAEIT'),
                      feat_reduce = './imputation/TEA_s1/scVAEIT/gene.csv')
B5_rna_df_2$task <- c('scVAEIT')

B3_rna_df <- rbind(B3_rna_df_1,B3_rna_df_2)
B5_rna_df <- rbind(B5_rna_df_1,B5_rna_df_2)
rna_df <- rbind(B3_rna_df,B5_rna_df)
rna_df

B3_adt_df <- evalue_adt(path_real = './PBMC_TEA/B3/',
                       path_impute = './imputation/TEA_s1',# 预测数据的路径
                       meta_path = './PBMC_TEA/B3/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'adt.rds',# 真实数据的名字
                       name_impute = c('/Palette/B3_adt.csv',
                                       '/MIDAS/B3_adt.csv',
                                       '/Multigrate/B3_adt.csv',
                                      '/scVAEIT/B3ADT.csv'),# 各个方法的分支路径
                       name_method = c('Palette','MIDAS','Multigrate','scVAEIT'))

B4_adt_df <- evalue_adt(path_real = './PBMC_TEA/B4/',
                       path_impute = './imputation/TEA_s1',# 预测数据的路径
                       meta_path = './PBMC_TEA/B4/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'adt.rds',# 真实数据的名字
                       name_impute = c('/Palette/B4_adt.csv',
                                       '/MIDAS/B4_adt.csv',
                                       '/Multigrate/B4_adt.csv',
                                      '/scVAEIT/B4ADT.csv'),# 各个方法的分支路径
                       name_method = c('Palette','MIDAS','Multigrate','scVAEIT'))

B3_adt_df$task <- c('Palette','MIDAS','Multigrate','scVAEIT')
B4_adt_df$task <- c('Palette','MIDAS','Multigrate','scVAEIT')

adt_df <- rbind(B3_adt_df,B4_adt_df)
adt_df

B4_atac_df_1 <- evalue_atac(path_real = './imputation/input/TEA/B4/',
                       path_impute = './imputation/TEA_s1',# 预测数据的路径
                       meta_path = './PBMC_TEA/B4/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'atac.rds',# 真实数据的名字
                       name_impute = c('Palette/B4_atac.mtx',
                                       'MIDAS/B4_atac.csv',
                                      'Multigrate/B4_atac.mtx'),# 各个方法的分支路径
                       name_method = c('Palette','MIDAS','Multigrate'))

B4_atac_df_2 <- evalue_atac(path_real = './imputation/input/TEA/B4/',
                       path_impute = './imputation/TEA_s1',# 预测数据的路径
                       meta_path = './PBMC_TEA/B4/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'atac.rds',# 真实数据的名字
                       name_impute = c('scVAEIT/B4ATAC.csv'),# 各个方法的分支路径
                       name_method = c('scVAEIT'),
                         feat_reduce = './imputation/TEA_s1/scVAEIT/peak.csv')
B4_atac_df_1$task <- c('Palette','MIDAS','Multigrate')
B4_atac_df_2$task <- c('scVAEIT')
B4_atac_df <- rbind(B4_atac_df_1,B4_atac_df_2)

B5_atac_df_1 <- evalue_atac(path_real = './imputation/input/TEA/B5/',
                       path_impute = './imputation/TEA_s1',# 预测数据的路径
                       meta_path = './PBMC_TEA/B5/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'atac.rds',# 真实数据的名字
                       name_impute = c('Palette/B5_atac.mtx',
                                       'MIDAS/B5_atac.csv',
                                      'Multigrate/B5_atac.mtx'),# 各个方法的分支路径
                       name_method = c('Palette','MIDAS','Multigrate'))

B5_atac_df_2 <- evalue_atac(path_real = './imputation/input/TEA/B5/',
                       path_impute = './imputation/TEA_s1',# 预测数据的路径
                       meta_path = './PBMC_TEA/B5/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'atac.rds',# 真实数据的名字
                       name_impute = c('scVAEIT/B5ATAC.csv'),# 各个方法的分支路径
                       name_method = c('scVAEIT'),
                         feat_reduce = './imputation/TEA_s1/scVAEIT/peak.csv')
B5_atac_df_1$task <- c('Palette','MIDAS','Multigrate')
B5_atac_df_2$task <- c('scVAEIT')
B5_atac_df <- rbind(B5_atac_df_1,B5_atac_df_2)

atac_df <- rbind(B4_atac_df,B5_atac_df)
atac_df

res <- list(RNA = rna_df,
           ADT = adt_df,
           ATAC = atac_df)
saveRDS(res,'./imputation/TEA_s1/metrics.rds')






