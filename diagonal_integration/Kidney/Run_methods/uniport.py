import numpy as np
import pandas as pd
from scipy.io import mmread
import matplotlib.pyplot as plt
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

atac = np.array(pyreadr.read_r('./ATAC/co.rds')[None]).T
meta_atac = pyreadr.read_r('./meta.rds')[None]
atac_adata = ad.AnnData(
    atac, dtype=np.float32,obs = meta_atac
)
rna = np.array(pyreadr.read_r('./RNA/X.rds')[None]).T
meta_rna = pyreadr.read_r('./RNA/meta.rds')[None]
rna_adata = ad.AnnData(
    rna, dtype=np.float32,obs=meta_rna
)

rna_adata.obs['domain_id'] = 0
rna_adata.obs['domain_id'] = rna_adata.obs['domain_id'].astype('category')
rna_adata.obs['source'] = 'RNA'
rna_adata.obs['source'] = rna_adata.obs['source'].astype('category')

atac_adata.obs['domain_id'] = 1
atac_adata.obs['domain_id'] = atac_adata.obs['domain_id'].astype('category')
atac_adata.obs['source'] = 'ATAC'
atac_adata.obs['source'] = atac_adata.obs['source'].astype('category')

adata_cm = rna_adata.concatenate(atac_adata, join='inner', batch_key='domain_id')

sc.pp.normalize_total(adata_cm)
sc.pp.log1p(adata_cm)
sc.pp.highly_variable_genes(adata_cm, n_top_genes=1565, inplace=False, subset=True)

up.batch_scale(adata_cm)

sc.pp.normalize_total(rna_adata)
sc.pp.log1p(rna_adata)
sc.pp.highly_variable_genes(rna_adata, n_top_genes=2000, inplace=False, subset=True)
up.batch_scale(rna_adata)

sc.pp.normalize_total(atac_adata)
sc.pp.log1p(atac_adata)
sc.pp.highly_variable_genes(atac_adata, n_top_genes=1565, inplace=False, subset=True)
up.batch_scale(atac_adata)

adata = up.Run(adatas=[rna_adata,atac_adata], adata_cm=adata_cm, lambda_s=1.0)

sc.pp.neighbors(adata, use_rep='latent')
sc.tl.umap(adata)
sc.pl.umap(adata, color=['modality'])
sc.pl.umap(adata, color=['label'])

latent = adata.obsm['latent']

if isinstance(latent, np.ndarray):
    latent = latent[:, :20]
    latent_df = pd.DataFrame(latent)
else:
    latent_df = latent.iloc[:, :20]


latent_df.to_csv('./kidney/uniport.csv', index=False)
