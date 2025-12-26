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

adt = pd.read_csv('./CODEX/codex.csv')
meta_adt = pyreadr.read_r('./CODEX/meta.rds')[None]
adt_adata = ad.AnnData(
    adt.to_numpy().T, dtype=np.float32,obs = meta_adt
)

rna = np.array(pyreadr.read_r('./RNA/rna.rds')[None]).T
meta_rna = pyreadr.read_r('./RNA/meta.rds')[None]
rna_adata = ad.AnnData(
    rna, dtype=np.float32,obs=meta_rna
)

rna_co = np.array(pyreadr.read_r('./RNA/co_all.rds')[None]).T
adt_co = pd.read_csv('./CODEX/co_all.csv')

rna_s_adata = ad.AnnData(
    rna_co, dtype=np.float32,obs = meta_rna
)

adt_s_adata = ad.AnnData(
    adt_co.to_numpy().T, dtype=np.float32,obs = meta_adt
)

del adt, meta_adt, rna, meta_rna, rna_co, adt_co

# process rna_shared
sc.pp.normalize_total(rna_s_adata)
sc.pp.log1p(rna_s_adata)
sc.pp.scale(rna_s_adata)

rna_s_adata = rna_s_adata.X.copy()
adt_s_adata = adt_s_adata.X.copy()

sc.pp.scale(rna_adata)

rna_active = rna_adata.X
protein_active = adt_adata.X

fusor = mf.model.Fusor(
    shared_arr1=adt_s_adata,
    shared_arr2=rna_s_adata,
    active_arr1=protein_active,
    active_arr2=rna_active,
    labels1=None,
    labels2=None
)

fusor.split_into_batches(max_outward_size=protein_active.shape[0],
                         matching_ratio=1,
                         metacell_size=1)

fusor.construct_graphs(
    n_neighbors1= 15,
    n_neighbors2= 15,
    svd_components1=20,
    svd_components2=20,
    verbose=True
)

fusor.find_initial_pivots(wt1=0.7, wt2=0.7,
    svd_components1=20, svd_components2=20)

fusor.refine_pivots(
    svd_components1=20, 
    svd_components2=20,
    cca_components=20,
    n_iters=3,
    randomized_svd=False, 
    svd_runs=1,
    verbose=True
)

fusor.filter_bad_matches()

protein_cca,rna_cca = fusor.get_embedding(
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

