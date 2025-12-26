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

s3_rna_df = cal_metrics(
    path ='./imputation/Retina/',
    data_path = ['Palette/LGS3OS_rna.csv','MIDAS/LGS3OS_rna.csv',
                  'Multigrate/LGS3OS_rna.csv','scVAEIT/LGS3OS_rna.csv'],              # list[str]：每个方法的数据文件相对 path 的路径
    method_name = ['Palette','MIDAS','Multigrate','scVAEIT'],            # list[str]：方法名（与 data_path 对齐）
    meta_path = './LGS3OS/meta.rds',              # RDS/CSV：包含真实标签
    label_key = 'cell_type__custom',              # 真实标签列名
    norm = [None,'log','log',None],                   # list[str|None]：每个方法的标准化方法
    dim_method = 'pca',             # 'pca' | 'lsi'
    dim = 50,                    # int：降维维度
    is_ATAC=False,
    clust_method='louvain'
)

d1_rna_df = cal_metrics(
    path ='./imputation/Retina/',
    data_path = ['Palette/LVG1OD_rna.csv','MIDAS/LVG1OD_rna.csv',
                  'Multigrate/LVG1OD_rna.csv','scVAEIT/LVG1OD_rna.csv'],              # list[str]：每个方法的数据文件相对 path 的路径
    method_name = ['Palette','MIDAS','Multigrate','scVAEIT'],            # list[str]：方法名（与 data_path 对齐）
    meta_path = './LVG1OD/meta.rds',              # RDS/CSV：包含真实标签
    label_key = 'cell_type__custom',              # 真实标签列名
    norm = [None,'log','log',None],                   # list[str|None]：每个方法的标准化方法
    dim_method = 'pca',             # 'pca' | 'lsi'
    dim = 50,                    # int：降维维度
    is_ATAC=False,
    clust_method='louvain'
)

s1_rna_df = cal_metrics(
    path ='./imputation/Retina/',
    data_path = ['Palette/LVG1OS_rna.csv','MIDAS/LVG1OS_rna.csv',
                  'Multigrate/LVG1OS_rna.csv','scVAEIT/LVG1OS_rna.csv'],              # list[str]：每个方法的数据文件相对 path 的路径
    method_name = ['Palette','MIDAS','Multigrate','scVAEIT'],            # list[str]：方法名（与 data_path 对齐）
    meta_path = './LVG1OS/meta.rds',              # RDS/CSV：包含真实标签
    label_key = 'cell_type__custom',              # 真实标签列名
    norm = [None,'log','log',None],                   # list[str|None]：每个方法的标准化方法
    dim_method = 'pca',             # 'pca' | 'lsi'
    dim = 50,                    # int：降维维度
    is_ATAC=False,
    clust_method='louvain'
)

rna_df = pd.concat([s3_rna_df,d1_rna_df,s1_rna_df], axis=0)
rna_df['Task'] = 'RNA'

d2_atac_df = cal_metrics(
    path = './imputation/Retina/',
    data_path = ['Palette/LGS2OD_atac.mtx','MIDAS/LGS2OD_atac.csv',
                  'Multigrate/LGS2OD_atac.mtx','scVAEIT/LGS2OD_atac.csv'],              # list[str]：每个方法的数据文件相对 path 的路径
    method_name = ['Palette','MIDAS','Multigrate','scVAEIT'],            # list[str]：方法名（与 data_path 对齐）
    meta_path = './LGS2OD/meta.rds',              # RDS/CSV：包含真实标签
    label_key = 'cell_type__custom',              # 真实标签列名
    norm = [None,'tfidf',None,'tfidf'],                   # list[str|None]：每个方法的标准化方法
    dim_method = 'lsi',             # 'pca' | 'lsi'
    dim = 50,                    # int：降维维度
    is_ATAC=True,
    clust_method='louvain'
)

s2_atac_df = cal_metrics(
    path = './imputation/Retina/',
    data_path = ['Palette/LGS2OS_atac.mtx','MIDAS/LGS2OS_atac.csv',
                  'Multigrate/LGS2OS_atac.mtx','scVAEIT/LGS2OS_atac.csv'],              # list[str]：每个方法的数据文件相对 path 的路径
    method_name = ['Palette','MIDAS','Multigrate','scVAEIT'],            # list[str]：方法名（与 data_path 对齐）
    meta_path = './LGS2OS/meta.rds',              # RDS/CSV：包含真实标签
    label_key = 'cell_type__custom',              # 真实标签列名
    norm = [None,'tfidf',None,'tfidf'],                   # list[str|None]：每个方法的标准化方法
    dim_method = 'lsi',             # 'pca' | 'lsi'
    dim = 50,                    # int：降维维度
    is_ATAC=True,
    clust_method='louvain'
)

d3_atac_df = cal_metrics(
    path = './imputation/Retina/',
    data_path = ['Palette/LGS3OD_atac.mtx','MIDAS/LGS3OD_atac.csv',
                  'Multigrate/LGS3OD_atac.mtx','scVAEIT/LGS3OD_atac.csv'],              # list[str]：每个方法的数据文件相对 path 的路径
    method_name = ['Palette','MIDAS','Multigrate','scVAEIT'],            # list[str]：方法名（与 data_path 对齐）
    meta_path = './LGS3OD/meta.rds',              # RDS/CSV：包含真实标签
    label_key = 'cell_type__custom',              # 真实标签列名
    norm = [None,'tfidf',None,'tfidf'],                   # list[str|None]：每个方法的标准化方法
    dim_method = 'lsi',             # 'pca' | 'lsi'
    dim = 50,                    # int：降维维度
    is_ATAC=True,
    clust_method='louvain'
)

atac_df = pd.concat([d2_atac_df,s2_atac_df,d3_atac_df], axis=0)
atac_df['Task'] = 'ATAC'

df = pd.concat([rna_df,atac_df], axis=0)
df.to_csv('./imputation/Retina/metrics_cluster.csv', index=False)
