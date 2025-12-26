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

merf = np.array(pyreadr.read_r('./MERFISH_ATAC/MERFISH/merfish.rds')[None]).T
meta_merf = pyreadr.read_r('./MERFISH_ATAC/MERFISH/meta.rds')[None]
merf_adata = ad.AnnData(
    merf, dtype=np.float32,obs = meta_merf
)

atac = np.array(pyreadr.read_r('./MERFISH_ATAC/ATAC/LSI.rds')[None])
meta_atac = pyreadr.read_r('./MERFISH_ATAC/ATAC/meta.rds')[None]
atac_adata = ad.AnnData(
    atac[:,0:30], dtype=np.float32,obs=meta_atac
)

merf_co = np.array(pyreadr.read_r('./MERFISH_ATAC/MERFISH/co100_RF.rds')[None]).T
atac_co = np.array(pyreadr.read_r('./MERFISH_ATAC/ATAC/co100_RF.rds')[None]).T

merf_s_adata = ad.AnnData(
    merf_co, dtype=np.float32,obs = meta_merf
)

atac_s_adata = ad.AnnData(
    atac_co, dtype=np.float32,obs = meta_atac
)
sc.pp.scale(merf_s_adata)
sc.pp.scale(atac_s_adata)
merf_s = merf_s_adata.X.copy()
atac_s = atac_s_adata.X.copy()
sc.pp.scale(merf_adata)

merf_active = merf_adata.X
atac_active = atac_adata.X

fusor = mf.model.Fusor(
    shared_arr1=atac_s,
    shared_arr2=merf_s,
    active_arr1=atac_active,
    active_arr2=merf_active,
    labels1=None,
    labels2=None
)

fusor.split_into_batches(max_outward_size=atac_active.shape[0],
                         matching_ratio=1,
                         metacell_size=1)

fusor.construct_graphs(
    n_neighbors1= 15,
    n_neighbors2= 15,
    svd_components1=20,
    svd_components2=20,
    verbose=True)

fusor.find_initial_pivots(svd_components1=20, svd_components2=20)

fusor.refine_pivots(
    svd_components1=20, 
    svd_components2=20,
    cca_components=20,
    n_iters=1,
    randomized_svd=False, 
    svd_runs=1,
    verbose=True
)

fusor.filter_bad_matches()

atac_cca, merf_cca = fusor.get_embedding(
    active_arr1=fusor.active_arr1,
    active_arr2=fusor.active_arr2)

cca_adata = ad.AnnData(
    np.concatenate((merf_cca, atac_cca), axis=0), 
    dtype=np.float32
)

latent_data = cca_adata.X
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    latent_df = latent_data

latent_df.to_csv('./Reduced_feat/RF100/maxfuse.csv', index=False)