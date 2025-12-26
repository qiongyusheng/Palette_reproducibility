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
sys.path.append('./metrics/')
from imputation_helper import cal_metrics

c1_rna_df = cal_metrics(
    path ='/imputation/BMMC_s/',
    data_path = ['Palette/cite1rna.csv','MIDAS/cite1rna.csv',
                 'Multigrate/cite1rna.csv','scVAEIT/cite1rna.csv'],              # list[str]：每个方法的数据文件相对 path 的路径
    method_name = ['Palette','MIDAS','Multigrate','scVAEIT'],            # list[str]：方法名（与 data_path 对齐）
    meta_path = '/BMMC_CITE_Multiome/CITE/s1d1_cite/meta.rds',              # RDS/CSV：包含真实标签
    label_key = 'celltype.l2',              # 真实标签列名
    norm = [None,'log','log',None],                   # list[str|None]：每个方法的标准化方法
    dim_method = 'pca',             # 'pca' | 'lsi'
    dim = 50,                    # int：降维维度
    is_ATAC=False,
    clust_method='louvain'
)

c2_rna_df = cal_metrics(
    path ='/imputation/BMMC_s/',
    data_path = ['Palette/cite2rna.csv','MIDAS/cite2rna.csv',
                 'Multigrate/cite2rna.csv','scVAEIT/cite2rna.csv'],              # list[str]：每个方法的数据文件相对 path 的路径
    method_name = ['Palette','MIDAS','Multigrate','scVAEIT'],            # list[str]：方法名（与 data_path 对齐）
    meta_path = '/BMMC_CITE_Multiome/CITE/s1d2_cite/meta.rds',              # RDS/CSV：包含真实标签
    label_key = 'celltype.l2',              # 真实标签列名
    norm = [None,'log','log',None],                   # list[str|None]：每个方法的标准化方法
    dim_method = 'pca',             # 'pca' | 'lsi'
    dim = 50,                    # int：降维维度
    is_ATAC=False,
    clust_method='louvain'
)

m2_rna_df = cal_metrics(
    path ='/imputation/BMMC_s/',
    data_path = ['Palette/mul2rna.csv','MIDAS/mul2rna.csv',
                 'Multigrate/mul2rna.csv','scVAEIT/mul2rna.csv'],              # list[str]：每个方法的数据文件相对 path 的路径
    method_name = ['Palette','MIDAS','Multigrate','scVAEIT'],            # list[str]：方法名（与 data_path 对齐）
    meta_path = '/BMMC_CITE_Multiome/Multiome/s1d2_multi/meta.rds',              # RDS/CSV：包含真实标签
    label_key = 'celltype.l2',              # 真实标签列名
    norm = [None,'log','log',None],                   # list[str|None]：每个方法的标准化方法
    dim_method = 'pca',             # 'pca' | 'lsi'
    dim = 50,                    # int：降维维度
    is_ATAC=False,
    clust_method='louvain'
)

m3_rna_df = cal_metrics(
    path ='/imputation/BMMC_s/',
    data_path = ['Palette/mul3rna.csv','MIDAS/mul3rna.csv',
                 'Multigrate/mul3rna.csv','scVAEIT/mul3rna.csv'],              # list[str]：每个方法的数据文件相对 path 的路径
    method_name = ['Palette','MIDAS','Multigrate','scVAEIT'],            # list[str]：方法名（与 data_path 对齐）
    meta_path = '/BMMC_CITE_Multiome/Multiome/s1d3_multi/meta.rds',              # RDS/CSV：包含真实标签
    label_key = 'celltype.l2',              # 真实标签列名
    norm = [None,'log','log',None],                   # list[str|None]：每个方法的标准化方法
    dim_method = 'pca',             # 'pca' | 'lsi'
    dim = 50,                    # int：降维维度
    is_ATAC=False,
    clust_method='louvain'
)

rna_df = pd.concat([c1_rna_df,c2_rna_df,m2_rna_df,m3_rna_df], axis=0)
rna_df['Task'] = 'RNA'

c1_adt_df = cal_metrics(
    path ='/imputation/BMMC_s/',
    data_path = ['Palette/cite1adt.csv','MIDAS/cite1adt.csv',
                 'Multigrate/cite1adt.csv','scVAEIT/cite1adt.csv'],              # list[str]：每个方法的数据文件相对 path 的路径
    method_name = ['Palette','MIDAS','Multigrate','scVAEIT'],            # list[str]：方法名（与 data_path 对齐）
    meta_path = '/BMMC_CITE_Multiome/CITE/s1d1_cite/meta.rds',              # RDS/CSV：包含真实标签
    label_key = 'celltype.l2',              # 真实标签列名
    norm = [None,'clr',None,None],                   # list[str|None]：每个方法的标准化方法
    dim_method = 'pca',             # 'pca' | 'lsi'
    dim = 50,                    # int：降维维度
    is_ATAC=False,
    clust_method='louvain'
)

c2_adt_df = cal_metrics(
    path ='/imputation/BMMC_s/',
    data_path = ['Palette/cite2adt.csv','MIDAS/cite2adt.csv',
                 'Multigrate/cite2adt.csv','scVAEIT/cite2adt.csv'],              # list[str]：每个方法的数据文件相对 path 的路径
    method_name = ['Palette','MIDAS','Multigrate','scVAEIT'],            # list[str]：方法名（与 data_path 对齐）
    meta_path = '/BMMC_CITE_Multiome/CITE/s1d2_cite/meta.rds',              # RDS/CSV：包含真实标签
    label_key = 'celltype.l2',              # 真实标签列名
    norm = [None,'clr',None,None],                   # list[str|None]：每个方法的标准化方法
    dim_method = 'pca',             # 'pca' | 'lsi'
    dim = 50,                    # int：降维维度
    is_ATAC=False,
    clust_method='louvain'
)

adt_df = pd.concat([c1_adt_df,c2_adt_df], axis=0)
adt_df['Task'] = 'ADT'

m2_atac_df = cal_metrics(
    path ='/imputation/BMMC_s/',
    data_path = ['Palette/mul2atac.mtx','MIDAS/mul2atac.csv',
                 'Multigrate/mul2atac.mtx','scVAEIT/mul2atac.csv'],              # list[str]：每个方法的数据文件相对 path 的路径
    method_name = ['Palette','MIDAS','Multigrate','scVAEIT'],            # list[str]：方法名（与 data_path 对齐）
    meta_path = '/imputation/input/BMMC_s/Multiome2atac/meta.rds',              # RDS/CSV：包含真实标签
    label_key = 'celltype.l2',              # 真实标签列名
    norm = [None,'tfidf',None,'tfidf'],                   # list[str|None]：每个方法的标准化方法
    dim_method = 'lsi',             # 'pca' | 'lsi'
    dim = 50,                    # int：降维维度
    is_ATAC=True,
    clust_method='louvain'
)

m3_atac_df = cal_metrics(
    path ='/imputation/BMMC_s/',
    data_path = ['Palette/mul3atac.mtx','MIDAS/mul3atac.csv',
                 'Multigrate/mul3atac.mtx','scVAEIT/mul3atac.csv'],              # list[str]：每个方法的数据文件相对 path 的路径
    method_name = ['Palette','MIDAS','Multigrate','scVAEIT'],            # list[str]：方法名（与 data_path 对齐）
    meta_path = '/imputation/input/BMMC_s/Multiome3atac/meta.rds',              # RDS/CSV：包含真实标签
    label_key = 'celltype.l2',              # 真实标签列名
    norm = [None,'tfidf',None,'tfidf'],                   # list[str|None]：每个方法的标准化方法
    dim_method = 'lsi',             # 'pca' | 'lsi'
    dim = 50,                    # int：降维维度
    is_ATAC=True,
    clust_method='louvain'
)

atac_df = pd.concat([m2_atac_df,m3_atac_df], axis=0)
atac_df['Task'] = 'ATAC'

df = pd.concat([rna_df,adt_df,atac_df], axis=0)

df.to_csv('/imputation/BMMC_s/metrics_cluster.csv', index=False)
