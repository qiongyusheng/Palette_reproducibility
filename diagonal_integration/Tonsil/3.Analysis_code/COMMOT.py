import os
import gc
import ot
import pickle
import anndata
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.io
from scipy import sparse
from scipy.stats import spearmanr, pearsonr
from scipy.spatial import distance_matrix
from scipy.sparse import coo_matrix
import matplotlib.pyplot as plt
import pyreadr
import commot as ct

data = scipy.io.mmread('./infer_rna=.mtx')
gene = pd.read_csv('./gene.csv')
bcd = pd.read_csv('./barcode.csv')
meta = pyreadr.read_r("./meta.rds")

adata = anndata.AnnData(X=data.T.tocsr(), obs=meta[None])

adata.var.index = gene['x']
adata.obs.index = bcd['x']

adata.obsm['spatial'] = meta[None][['X', 'Y']].to_numpy()
adata.var['gene_ids'] = gene['x']
df_cellphonedb = ct.pp.ligand_receptor_database(species='human',
                                             signaling_type='Secreted Signaling', 
                                             database='CellPhoneDB_v4.0')

unique_0 = df_cellphonedb.iloc[:, 0].unique()
unique_1 = df_cellphonedb.iloc[:, 1].unique()

final_unique = np.unique(np.concatenate([unique_0, unique_1]))

intersection = np.intersect1d(final_unique, gene['x'])

adata_sub = adata[:,intersection]
del data
del adata
import gc
gc.collect()

df_cellphonedb_filtered = ct.pp.filter_lr_database(df_cellphonedb, adata_sub, min_cell_pct=0.01)

ct.tl.spatial_communication(adata_sub,
                            database_name='CellPhoneDB_v4.0', 
                            df_ligrec=df_cellphonedb_filtered, 
                            dis_thr=500, 
                            heteromeric=True, 
                            pathway_sum=True)

adata_sub.var['gene_ids'] = adata_sub.var['gene_ids'].astype(str)
adata_sub.write('./result.h5ad')