library(tools)
library(glue)
library(data.table)
library(Seurat)
library(gridExtra)
library(SingleCellExperiment)
library(parallel)
library(scran)
source('./metrics/imputation_helper.R')

c1_rna_df_1 <- evalue_rna(path_real = './BMMC_CITE_Multiome/CITE/s1d1_cite',
                       path_impute = './imputation/BMMC_s',# 预测数据的路径
                       meta_path = './BMMC_CITE_Multiome/CITE/s1d1_cite/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'rna.rds',# 真实数据的名字
                       name_impute = c('Palette/cite1rna.csv',
                                       'MIDAS/cite1rna.csv',
                                       'Multigrate/cite1rna.csv'),# 各个方法的分支路径
                       name_method = c('Palette','MIDAS','Multigrate'))
c1_rna_df_1$task <- c('Palette','MIDAS','Multigrate')

c1_rna_df_2 <- evalue_rna(path_real = './BMMC_CITE_Multiome/CITE/s1d1_cite',
                       path_impute = './imputation/BMMC_s',# 预测数据的路径
                       meta_path = './BMMC_CITE_Multiome/CITE/s1d1_cite/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'rna.rds',# 真实数据的名字
                       name_impute = c('/scVAEIT/cite1rna.csv'),# 各个方法的分支路径
                       name_method = c('scVAEIT'),
                      feat_reduce = './imputation/BMMC_s/scVAEIT/gene.csv')
c1_rna_df_2$task <- c('scVAEIT')

c1_rna_df <- rbind(c1_rna_df_1,c1_rna_df_2)

c2_rna_df_1 <- evalue_rna(path_real = './BMMC_CITE_Multiome/CITE/s1d2_cite',
                       path_impute = './imputation/BMMC_s',# 预测数据的路径
                       meta_path = './BMMC_CITE_Multiome/CITE/s1d2_cite/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'rna.rds',# 真实数据的名字
                       name_impute = c('Palette/cite2rna.csv',
                                       'MIDAS/cite2rna.csv',
                                       'Multigrate/cite2rna.csv'),# 各个方法的分支路径
                       name_method = c('Palette','MIDAS','Multigrate'))
c2_rna_df_1$task <- c('Palette','MIDAS','Multigrate')

c2_rna_df_2 <- evalue_rna(path_real = './BMMC_CITE_Multiome/CITE/s1d2_cite',
                       path_impute = './imputation/BMMC_s',# 预测数据的路径
                       meta_path = './BMMC_CITE_Multiome/CITE/s1d2_cite/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'rna.rds',# 真实数据的名字
                       name_impute = c('/scVAEIT/cite2rna.csv'),# 各个方法的分支路径
                       name_method = c('scVAEIT'),
                      feat_reduce = './imputation/BMMC_s/scVAEIT/gene.csv')
c2_rna_df_2$task <- c('scVAEIT')

c2_rna_df <- rbind(c2_rna_df_1,c2_rna_df_2)

m2_rna_df_1 <- evalue_rna(path_real = './BMMC_CITE_Multiome/Multiome/s1d2_multi',
                       path_impute = './imputation/BMMC_s',# 预测数据的路径
                       meta_path = './BMMC_CITE_Multiome/Multiome/s1d2_multi/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'rna.rds',# 真实数据的名字
                       name_impute = c('/Palette/mul2rna.csv',
                                       '/MIDAS/mul2rna.csv',
                                       '/Multigrate/mul2rna.csv'),# 各个方法的分支路径
                       name_method = c('Palette','MIDAS','Multigrate'))
m2_rna_df_1$task <- c('Palette','MIDAS','Multigrate')

m2_rna_df_2 <- evalue_rna(path_real = './BMMC_CITE_Multiome/Multiome/s1d2_multi',
                       path_impute = './imputation/BMMC_s',# 预测数据的路径
                       meta_path = './BMMC_CITE_Multiome/Multiome/s1d2_multi/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'rna.rds',# 真实数据的名字
                       name_impute = c('/scVAEIT/mul2rna.csv'),# 各个方法的分支路径
                       name_method = c('scVAEIT'),
                      feat_reduce = './imputation/BMMC_s/scVAEIT/gene.csv')
m2_rna_df_2$task <- c('scVAEIT')

m2_rna_df <- rbind(m2_rna_df_1,m2_rna_df_2)

m3_rna_df_1 <- evalue_rna(path_real = './BMMC_CITE_Multiome/Multiome/s1d3_multi',
                       path_impute = './imputation/BMMC_s',# 预测数据的路径
                       meta_path = './BMMC_CITE_Multiome/Multiome/s1d3_multi/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'rna.rds',# 真实数据的名字
                       name_impute = c('/Palette/mul3rna.csv',
                                       '/MIDAS/mul3rna.csv',
                                       '/Multigrate/mul3rna.csv'),# 各个方法的分支路径
                       name_method = c('Palette','MIDAS','Multigrate'))
m3_rna_df_1$task <- c('Palette','MIDAS','Multigrate')

m3_rna_df_2 <- evalue_rna(path_real = './BMMC_CITE_Multiome/Multiome/s1d3_multi',
                       path_impute = './imputation/BMMC_s',# 预测数据的路径
                       meta_path = './BMMC_CITE_Multiome/Multiome/s1d3_multi/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'rna.rds',# 真实数据的名字
                       name_impute = c('/scVAEIT/mul3rna.csv'),# 各个方法的分支路径
                       name_method = c('scVAEIT'),
                      feat_reduce = './imputation/BMMC_s/scVAEIT/gene.csv')
m3_rna_df_2$task <- c('scVAEIT')

m3_rna_df <- rbind(m3_rna_df_1,m3_rna_df_2)

rna_df <- rbind(c1_rna_df,c2_rna_df,m2_rna_df,m3_rna_df)
rna_df

c1_adt_df <- evalue_adt(path_real = './BMMC_CITE_Multiome/CITE/s1d1_cite',
                       path_impute = './imputation/BMMC_s',# 预测数据的路径
                       meta_path = './BMMC_CITE_Multiome/CITE/s1d1_cite/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'adt.rds',# 真实数据的名字
                       name_impute = c('/Palette/cite1adt.csv',
                                       '/MIDAS/cite1adt.csv',
                                       '/Multigrate/cite1adt.csv',
                                      '/scVAEIT/cite1adt.csv'),# 各个方法的分支路径
                       name_method = c('Palette','MIDAS','Multigrate','scVAEIT'))

c2_adt_df <- evalue_adt(path_real = './BMMC_CITE_Multiome/CITE/s1d2_cite',
                       path_impute = './imputation/BMMC_s',# 预测数据的路径
                       meta_path = './BMMC_CITE_Multiome/CITE/s1d2_cite/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'adt.rds',# 真实数据的名字
                       name_impute = c('/Palette/cite2adt.csv',
                                       '/MIDAS/cite2adt.csv',
                                       '/Multigrate/cite2adt.csv',
                                      '/scVAEIT/cite2adt.csv'),# 各个方法的分支路径
                       name_method = c('Palette','MIDAS','Multigrate','scVAEIT'))

c1_adt_df$task <- c('Palette','MIDAS','Multigrate','scVAEIT')
c2_adt_df$task <- c('Palette','MIDAS','Multigrate','scVAEIT')

adt_df <- rbind(c1_adt_df,c2_adt_df)
adt_df

m2_atac_df_1 <- evalue_atac(path_real = './imputation/input/BMMC_s/Multiome2atac',
                       path_impute = './imputation/BMMC_s',# 预测数据的路径
                       meta_path = './imputation/input/BMMC_s/Multiome2atac/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'atac.rds',# 真实数据的名字
                       name_impute = c('Palette/mul2atac.mtx',
                                       'MIDAS/mul2atac.csv'),# 各个方法的分支路径
                       name_method = c('Palette','MIDAS'))

m2_atac_df_2 <- evalue_atac(path_real = './imputation/input/BMMC_s/Multiome2atac',
                       path_impute = './imputation/BMMC_s',# 预测数据的路径
                       meta_path = './imputation/input/BMMC_s/Multiome2atac/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'atac.rds',# 真实数据的名字
                       name_impute = c('Multigrate/mul2atac.mtx'),# 各个方法的分支路径
                       name_method = c('Multigrate'))

m2_atac_df_3 <- evalue_atac(path_real = './imputation/input/BMMC_s/Multiome2atac',
                       path_impute = './imputation/BMMC_s',# 预测数据的路径
                       meta_path = './imputation/input/BMMC_s/Multiome2atac/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'atac.rds',# 真实数据的名字
                       name_impute = c('scVAEIT/mul2atac.csv'),# 各个方法的分支路径
                       name_method = c('scVAEIT'),
                       feat_reduce = './imputation/BMMC_s/scVAEIT/peak.csv')
m2_atac_df_1$task <- c('Palette','MIDAS')
m2_atac_df_2$task <- c('Multigrate')
m2_atac_df_3$task <- c('scVAEIT')
m2_atac_df <- rbind(m2_atac_df_1,m2_atac_df_2,m2_atac_df_3)

m3_atac_df_1 <- evalue_atac(path_real = './imputation/input/BMMC_s/Multiome3atac',
                       path_impute = './imputation/BMMC_s',# 预测数据的路径
                       meta_path = './imputation/input/BMMC_s/Multiome3atac/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'atac.rds',# 真实数据的名字
                       name_impute = c('Palette/mul3atac.mtx',
                                       'MIDAS/mul3atac.csv'),# 各个方法的分支路径
                       name_method = c('Palette','MIDAS'))

m3_atac_df_2 <- evalue_atac(path_real = './imputation/input/BMMC_s/Multiome3atac',
                       path_impute = './imputation/BMMC_s',# 预测数据的路径
                       meta_path = './imputation/input/BMMC_s/Multiome3atac/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'atac.rds',# 真实数据的名字
                       name_impute = c('Multigrate/mul3atac.mtx'),# 各个方法的分支路径
                       name_method = c('Multigrate'))


m3_atac_df_3 <- evalue_atac(path_real = './imputation/input/BMMC_s/Multiome3atac',
                       path_impute = './imputation/BMMC_s',# 预测数据的路径
                       meta_path = './imputation/input/BMMC_s/Multiome3atac/meta.rds',# meta数据路径
                       group = 'celltype.l2',# 使用的meta列名
                       name_real = 'atac.rds',# 真实数据的名字
                       name_impute = c('scVAEIT/mul3atac.csv'),# 各个方法的分支路径
                       name_method = c('scVAEIT'),
                       feat_reduce = './imputation/BMMC_s/scVAEIT/peak.csv')
m3_atac_df_1$task <- c('Palette','MIDAS')
m3_atac_df_2$task <- c('Multigrate')
m3_atac_df_3$task <- c('scVAEIT')
m3_atac_df <- rbind(m3_atac_df_1,m3_atac_df_2,m3_atac_df_3)

atac_df <- rbind(m2_atac_df,m3_atac_df)
atac_df

res <- list(RNA = rna_df,
           ADT = adt_df,
           ATAC = atac_df)
saveRDS(res,'./imputation/BMMC_s/metrics.rds')
