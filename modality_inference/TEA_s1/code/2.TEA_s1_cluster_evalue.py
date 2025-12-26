import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.io as sio
from scipy.sparse import csr_matrix, issparse
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.decomposition import TruncatedSVD
from sklearn import metrics
import muon
import episcanpy as epi

import sys
sys.path.append('./metrics')
from imputation_helper import cal_metrics

B3_rna_df = cal_metrics(
    path ='/imputation/TEA_s1/',
    data_path = ['Palette/B3_rna.csv','MIDAS/B3_rna.csv',
                 'Multigrate/B3_rna.csv','scVAEIT/B3RNA.csv'],              # list[str]：每个方法的数据文件相对 path 的路径
    method_name = ['Palette','MIDAS','Multigrate','scVAEIT'],            # list[str]：方法名（与 data_path 对齐）
    meta_path = './PBMC_TEA/B3/meta.rds',              # RDS/CSV：包含真实标签
    label_key = 'celltype.l2',              # 真实标签列名
    norm = [None,'log','log',None],                   # list[str|None]：每个方法的标准化方法
    dim_method = 'pca',             # 'pca' | 'lsi'
    dim = 50,                    # int：降维维度
    is_ATAC=False,
    clust_method='louvain'
)

B5_rna_df = cal_metrics(
    path ='/imputation/TEA_s1/',
    data_path = ['Palette/B5_rna.csv','MIDAS/B5_rna.csv',
                 'Multigrate/B5_rna.csv','scVAEIT/B5RNA.csv'],              # list[str]：每个方法的数据文件相对 path 的路径
    method_name = ['Palette','MIDAS','Multigrate','scVAEIT'],            # list[str]：方法名（与 data_path 对齐）
    meta_path = './PBMC_TEA/B5/meta.rds',              # RDS/CSV：包含真实标签
    label_key = 'celltype.l2',              # 真实标签列名
    norm = [None,'log','log',None],                   # list[str|None]：每个方法的标准化方法
    dim_method = 'pca',             # 'pca' | 'lsi'
    dim = 50,                    # int：降维维度
    is_ATAC=False,
    clust_method='louvain'
)
rna_df = pd.concat([B3_rna_df,B5_rna_df], axis=0)
rna_df['Task'] = 'RNA'

B3_adt_df = cal_metrics(
    path ='/imputation/TEA_s1/',
    data_path = ['Palette/B3_adt.csv','MIDAS/B3_adt.csv',
                 'Multigrate/B3_adt.csv','scVAEIT/B3ADT.csv'],              # list[str]：每个方法的数据文件相对 path 的路径
    method_name = ['Palette','MIDAS','Multigrate','scVAEIT'],            # list[str]：方法名（与 data_path 对齐）
    meta_path = './PBMC_TEA/B3/meta.rds',              # RDS/CSV：包含真实标签
    label_key = 'celltype.l2',              # 真实标签列名
    norm = [None,'clr',None,None],                   # list[str|None]：每个方法的标准化方法
    dim_method = 'pca',             # 'pca' | 'lsi'
    dim = 30,                    # int：降维维度
    is_ATAC=False,
    clust_method='louvain'
)

B4_adt_df = cal_metrics(
    path ='/imputation/TEA_s1/',
    data_path = ['Palette/B4_adt.csv','MIDAS/B4_adt.csv',
                 'Multigrate/B4_adt.csv','scVAEIT/B4ADT.csv'],              # list[str]：每个方法的数据文件相对 path 的路径
    method_name = ['Palette','MIDAS','Multigrate','scVAEIT'],            # list[str]：方法名（与 data_path 对齐）
    meta_path = './PBMC_TEA/B4/meta.rds',              # RDS/CSV：包含真实标签
    label_key = 'celltype.l2',              # 真实标签列名
    norm = [None,'clr',None,None],                   # list[str|None]：每个方法的标准化方法
    dim_method = 'pca',             # 'pca' | 'lsi'
    dim = 30,                    # int：降维维度
    is_ATAC=False,
    clust_method='louvain'
)
adt_df = pd.concat([B3_adt_df,B4_adt_df], axis=0)
adt_df['Task'] = 'ADT'

B4_atac_df = cal_metrics(
    path ='/imputation/TEA_s1/',
    data_path = ['Palette/B4_atac.mtx','MIDAS/B4_atac.csv',
                 'Multigrate/B4_atac.mtx','scVAEIT/B4ATAC.csv'],              # list[str]：每个方法的数据文件相对 path 的路径
    method_name = ['Palette','MIDAS','Multigrate','scVAEIT'],            # list[str]：方法名（与 data_path 对齐）
    meta_path = './PBMC_TEA/B4/meta.rds',              # RDS/CSV：包含真实标签
    label_key = 'celltype.l2',              # 真实标签列名
    norm = [None,'tfidf',None,'tfidf'],                 # list[str|None]：每个方法的标准化方法
    dim_method = 'lsi',             # 'pca' | 'lsi'
    dim = 50,                    # int：降维维度
    is_ATAC=True,
    clust_method='louvain'
)

B5_atac_df = cal_metrics(
    path ='/imputation/TEA_s1/',
    data_path = ['Palette/B5_atac.mtx','MIDAS/B5_atac.csv',
                 'Multigrate/B5_atac.mtx','scVAEIT/B5ATAC.csv'],              # list[str]：每个方法的数据文件相对 path 的路径
    method_name = ['Palette','MIDAS','Multigrate','scVAEIT'],            # list[str]：方法名（与 data_path 对齐）
    meta_path = './PBMC_TEA/B5/meta.rds',              # RDS/CSV：包含真实标签
    label_key = 'celltype.l2',              # 真实标签列名
    norm = [None,'tfidf',None,'tfidf'],                 # list[str|None]：每个方法的标准化方法
    dim_method = 'lsi',             # 'pca' | 'lsi'
    dim = 50,                    # int：降维维度
    is_ATAC=True,
    clust_method='louvain'
)

atac_df = pd.concat([B4_atac_df,B5_atac_df], axis=0)
atac_df['Task'] = 'ATAC'

df = pd.concat([rna_df,adt_df,atac_df], axis=0)
df.to_csv('/imputation/TEA_s1/metrics_cluster.csv', index=False)
