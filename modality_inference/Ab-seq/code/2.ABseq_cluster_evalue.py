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

a2_rna_df = cal_metrics(
    path = './ABseq/predict/',
    data_path = ['palette/Age2rna_k10.csv','midas/Age2rna.csv',
                  'multigrate/Age2rna.csv','scvaeit/Age2rna.csv'],              # list[str]：每个方法的数据文件相对 path 的路径
    method_name = ['Palette','MIDAS','Multigrate','scVAEIT'],            # list[str]：方法名（与 data_path 对齐）
    meta_path = './Age2/meta.rds',              # RDS/CSV：包含真实标签
    label_key = 'celltype.l2',              # 真实标签列名
    norm = [None,'log','log',None],                   # list[str|None]：每个方法的标准化方法
    dim_method = 'pca',             # 'pca' | 'lsi'
    dim = 50,                    # int：降维维度
    is_ATAC=False,
    clust_method='louvain'
)

m3_rna_df = cal_metrics(
    path = './ABseq/predict/',
    data_path = ['palette/BM3rna_k10.csv','midas/BM3rna.csv',
                  'multigrate/BM3rna.csv','scvaeit/BM3rna.csv'],              # list[str]：每个方法的数据文件相对 path 的路径
    method_name = ['Palette','MIDAS','Multigrate','scVAEIT'],            # list[str]：方法名（与 data_path 对齐）
    meta_path = './BM3/meta.rds',              # RDS/CSV：包含真实标签
    label_key = 'celltype.l2',              # 真实标签列名
    norm = [None,'log','log',None],                   # list[str|None]：每个方法的标准化方法
    dim_method = 'pca',             # 'pca' | 'lsi'
    dim = 50,                    # int：降维维度
    is_ATAC=False,
    clust_method='louvain'
)
rna_df = pd.concat([a2_rna_df,m3_rna_df], axis=0)
rna_df['Task'] = 'RNA'


a3_adt_df = cal_metrics(
    path = './ABseq/predict/',
    data_path = ['palette/Age3adt_k10.csv','midas/Age3adt.csv',
                  'multigrate/Age3adt.csv','scvaeit/Age3adt.csv'],              # list[str]：每个方法的数据文件相对 path 的路径
    method_name = ['Palette','MIDAS','Multigrate','scVAEIT'],            # list[str]：方法名（与 data_path 对齐）
    meta_path = './Age3/meta.rds',              # RDS/CSV：包含真实标签
    label_key = 'celltype.l2',              # 真实标签列名
    norm = [None,'clr',None,None],                   # list[str|None]：每个方法的标准化方法
    dim_method = 'pca',             # 'pca' | 'lsi'
    dim = 50,                    # int：降维维度
    is_ATAC=False,
    clust_method='louvain'
)

m2_adt_df = cal_metrics(
    path = './ABseq/predict/',
    data_path = ['palette/BM2adt_k10.csv','midas/BM2adt.csv',
                  'multigrate/BM2adt.csv','scvaeit/BM2adt.csv'],              # list[str]：每个方法的数据文件相对 path 的路径
    method_name = ['Palette','MIDAS','Multigrate','scVAEIT'],            # list[str]：方法名（与 data_path 对齐）
    meta_path = './BM2/meta.rds',              # RDS/CSV：包含真实标签
    label_key = 'celltype.l2',              # 真实标签列名
    norm = [None,'clr',None,None],                   # list[str|None]：每个方法的标准化方法
    dim_method = 'pca',             # 'pca' | 'lsi'
    dim = 50,                    # int：降维维度
    is_ATAC=False,
    clust_method='louvain'
)

adt_df = pd.concat([a3_adt_df,m2_adt_df], axis=0)
adt_df['Task'] = 'ADT'

df = pd.concat([rna_df,adt_df], axis=0)
df.to_csv('./ABseq/metrics_cluster.csv', index=False)