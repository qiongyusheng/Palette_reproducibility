import numpy as np
import pandas as pd
from scipy.io import mmread
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (6, 4)

from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

import anndata as ad
import scanpy as sc
import maxfuse as mf
import scanpy as sc
import anndata as ad
import json
import os
import pyreadr

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

# process rna_shared
sc.pp.normalize_total(rna_s_adata)
sc.pp.log1p(rna_s_adata)
sc.pp.scale(rna_s_adata)

rna_s_adata = rna_s_adata.X.copy()

# process protein_shared
sc.pp.normalize_total(adt_s_adata)
sc.pp.log1p(adt_s_adata)
sc.pp.scale(adt_s_adata)

adt_s_adata = adt_s_adata.X.copy()

# process all RNA features
sc.pp.normalize_total(rna_adata)
sc.pp.log1p(rna_adata)
#sc.pp.highly_variable_genes(rna_adata)
# only retain highly variable genes
#rna_adata = rna_adata[:, rna_adata.var.highly_variable].copy()
sc.pp.scale(rna_adata)

# process all protein features
sc.pp.normalize_total(adt_adata)
sc.pp.log1p(adt_adata)
sc.pp.scale(adt_adata)

rna_active = rna_adata.X
protein_active = adt_adata.X

# inspect shape of the four matrices
print(rna_active.shape)
print(protein_active.shape)
print(rna_s_adata.shape)
print(adt_s_adata.shape)

# call constructor for Fusor object
# which is the main object for running MaxFuse pipeline
fusor = mf.model.Fusor(
    shared_arr1=rna_s_adata,
    shared_arr2=adt_s_adata,
    active_arr1=rna_active,
    active_arr2=protein_active,
    labels1=None,
    labels2=None
)

fusor.split_into_batches()

fusor.construct_graphs(
    n_neighbors1= 15,
    n_neighbors2= 15,
    svd_components1=20,
    svd_components2=20,
    verbose=True
)
fusor.find_initial_pivots(svd_components1=20, svd_components2=20)

fusor.refine_pivots(svd_components1=20, svd_components2=None,
    cca_components=20,
    n_iters=1,
    verbose=True
)

fusor.filter_bad_matches(target='pivot', filter_prop=0)

rna_cca, protein_cca = fusor.get_embedding(
    active_arr1=fusor.active_arr1,
    active_arr2=fusor.active_arr2
)


cca_adata = ad.AnnData(
    np.concatenate((rna_cca, protein_cca), axis=0), 
    dtype=np.float32
)


latent_data = cca_adata.X
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    latent_df = latent_data

latent_df.to_csv('./maxfuse.csv', index=False)
