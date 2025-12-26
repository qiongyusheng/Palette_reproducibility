import numpy as np
import pandas as pd
from scipy.io import mmread
import matplotlib.pyplot as plt
from scipy import sparse
plt.rcParams["figure.figsize"] = (6, 4)

from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

import anndata as ad
import scanpy as sc
import anndata as ad
import json
import os
import pyreadr
import uniport as up
print(up.__version__)

rds_data = pyreadr.read_r('./MERFISH_ATAC/ATAC/rna.rds')
atac = np.array(rds_data[None]).T
atac_sparse = sparse.csr_matrix(atac)
del atac
gene_name = np.array(rds_data[None].index, dtype=str)
meta_atac = pyreadr.read_r('./MERFISH_ATAC/ATAC/meta.rds')[None]
atac_adata = ad.AnnData(
    atac_sparse, dtype=np.float32,obs = meta_atac
)
atac_adata.var_names = gene_name

rds_data = pyreadr.read_r('./MERFISH_ATAC/MERFISH/X.rds')
rna = np.array(rds_data[None]).T
gene_name = np.array(rds_data[None].index, dtype=str)
meta_rna = pyreadr.read_r('./MERFISH_ATAC/MERFISH/meta.rds')[None]
rna_adata = ad.AnnData(
    rna, dtype=np.float32,obs=meta_rna
)
rna_adata.var_names = gene_name

del rds_data
del gene_name
del meta_rna
del meta_atac
import gc
gc.collect()

rna_adata.obs['domain_id'] = 0
rna_adata.obs['domain_id'] = rna_adata.obs['domain_id'].astype('category')
rna_adata.obs['source'] = 'RNA'
rna_adata.obs['source'] = rna_adata.obs['source'].astype('category')

atac_adata.obs['domain_id'] = 1
atac_adata.obs['domain_id'] = atac_adata.obs['domain_id'].astype('category')
atac_adata.obs['source'] = 'ATAC'
atac_adata.obs['source'] = atac_adata.obs['source'].astype('category')

del rna
del atac_sparse
gc.collect()

adata_cm = rna_adata.concatenate(atac_adata, join='inner', batch_key='domain_id')

sc.pp.normalize_total(adata_cm)
sc.pp.log1p(adata_cm)
sc.pp.highly_variable_genes(adata_cm, n_top_genes=253, inplace=False, subset=True)
up.batch_scale(adata_cm)

sc.pp.normalize_total(rna_adata)
sc.pp.log1p(rna_adata)
sc.pp.highly_variable_genes(rna_adata, n_top_genes=254, inplace=False, subset=True)
up.batch_scale(rna_adata)

sc.pp.normalize_total(atac_adata)
sc.pp.log1p(atac_adata)
sc.pp.highly_variable_genes(atac_adata, n_top_genes=2000, inplace=False, subset=True)
up.batch_scale(atac_adata)

adata = up.Run(adatas=[rna_adata,atac_adata], adata_cm=adata_cm, lambda_s=1.0)

latent_data = adata.obsm['latent']
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    latent_df = latent_data

latent_df.to_csv('./uniport.csv', index=False)