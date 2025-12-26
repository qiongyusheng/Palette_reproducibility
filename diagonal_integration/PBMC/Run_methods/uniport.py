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

adt = np.array(pyreadr.read_r('./cytof/cytof.rds')[None]).T
meta_adt = pyreadr.read_r('./cytof/meta.rds')[None]
adt_adata = ad.AnnData(
    adt, dtype=np.float32,obs = meta_adt
)

rna = np.array(pyreadr.read_r('./RNA/rna.rds')[None]).T
meta_rna = pyreadr.read_r('./RNA/meta.rds')[None]
rna_adata = ad.AnnData(
    rna, dtype=np.float32,obs=meta_rna
)

rna_co = np.array(pyreadr.read_r('./RNA/co.rds')[None]).T
adt_co = np.array(pyreadr.read_r('./cytof/co.rds')[None]).T

rna_s_adata = ad.AnnData(
    rna_co, dtype=np.float32,obs = meta_rna
)

adt_s_adata = ad.AnnData(
    adt_co, dtype=np.float32,obs = meta_adt
)

rna_adata.obs['domain_id'] = 0
rna_adata.obs['domain_id'] = rna_adata.obs['domain_id'].astype('category')
rna_s_adata.obs['domain_id'] = 0
rna_s_adata.obs['domain_id'] = rna_s_adata.obs['domain_id'].astype('category')

rna_adata.obs['source'] = 'RNA'
rna_adata.obs['source'] = rna_adata.obs['source'].astype('category')
rna_s_adata.obs['source'] = 'RNA'
rna_s_adata.obs['source'] = rna_s_adata.obs['source'].astype('category')

adt_adata.obs['domain_id'] = 1
adt_adata.obs['domain_id'] = adt_adata.obs['domain_id'].astype('category')
adt_s_adata.obs['domain_id'] = 1
adt_s_adata.obs['domain_id'] = adt_s_adata.obs['domain_id'].astype('category')

adt_adata.obs['source'] = 'CyTOF'
adt_adata.obs['source'] = adt_adata.obs['source'].astype('category')
adt_s_adata.obs['source'] = 'CyTOF'
adt_s_adata.obs['source'] = adt_s_adata.obs['source'].astype('category')

adata_cm = rna_s_adata.concatenate(adt_s_adata, join='inner', batch_key='domain_id')

sc.pp.normalize_total(adata_cm)
sc.pp.log1p(adata_cm)
sc.pp.highly_variable_genes(adata_cm, n_top_genes=35, inplace=False, subset=True)

up.batch_scale(adata_cm)

sc.pp.normalize_total(rna_adata)
sc.pp.log1p(rna_adata)
sc.pp.highly_variable_genes(rna_adata, n_top_genes=2000, inplace=False, subset=True)
up.batch_scale(rna_adata)

sc.pp.normalize_total(adt_adata)
sc.pp.log1p(adt_adata)
sc.pp.highly_variable_genes(adt_adata, n_top_genes=30, inplace=False, subset=True)
up.batch_scale(adt_adata)

adata = up.Run(adatas=[rna_adata,adt_adata], adata_cm=adata_cm, lambda_s=1.0)


latent_data = adata.obsm['latent']

if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:

    latent_df = latent_data

latent_df.to_csv('./uniport.csv', index=False)

