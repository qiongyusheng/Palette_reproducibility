import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from scipy.spatial.distance import cdist
import muon as mu
import numpy as np
import anndata as ad
import scanpy as sc
import pyreadr
import scconfluence
import pandas as pd
from scipy.sparse import csr_matrix

merf = np.array(pyreadr.read_r('./MERFISH_ATAC/MERFISH/merfish.rds')[None]).T
meta_merf = pyreadr.read_r('./MERFISH_ATAC/MERFISH/meta.rds')[None]
merf_adata = ad.AnnData(
    merf, dtype=np.float32,obs = meta_merf
)

atac = np.array(pyreadr.read_r('./MERFISH_ATAC/ATAC/LSI.rds')[None])
meta_atac = pyreadr.read_r('./MERFISH_ATAC/ATAC/meta.rds')[None]
atac_adata = ad.AnnData(
    atac, dtype=np.float32,obs=meta_atac
)

merf_co = np.array(pyreadr.read_r('./MERFISH_ATAC/MERFISH/co150_RF.rds')[None]).T
atac_co = np.array(pyreadr.read_r('./MERFISH_ATAC/ATAC/co150_RF.rds')[None]).T

merf_s_adata = ad.AnnData(
    merf_co, dtype=np.float32,obs = meta_merf
)

atac_s_adata = ad.AnnData(
    atac_co, dtype=np.float32,obs = meta_atac
)

merf_adata.obs_names = [f"rna{i+1}" for i in range(merf_s_adata.n_obs)]
atac_adata.obs_names = [f"atac{i+1}" for i in range(atac_s_adata.n_obs)]

merf_s_adata.obs_names = [f"rna{i+1}" for i in range(merf_s_adata.n_obs)]
atac_s_adata.obs_names = [f"atac{i+1}" for i in range(atac_s_adata.n_obs)]

mdata = mu.MuData({"rna": merf_adata, "atac": atac_adata})

sc.pp.normalize_total(merf_s_adata, target_sum=10000.)
sc.pp.log1p(merf_s_adata)

# Log-normalize scATAC gene activities
sc.pp.normalize_total(atac_s_adata, target_sum=10000.)
sc.pp.log1p(atac_s_adata)

# Select highly variable genes for both modalities using the more reliable scRNA counts
cm_hvg_genes = sc.pp.highly_variable_genes(merf_s_adata, n_top_genes=1565, subset=False, inplace=False)
merf_s_adata = merf_s_adata[:, cm_hvg_genes["highly_variable"]].copy()
atac_s_adata = atac_s_adata[:, cm_hvg_genes["highly_variable"]].copy()

sc.pp.scale(merf_s_adata)
sc.pp.scale(atac_s_adata)

del merf,atac,merf_co,atac_co,merf_adata,atac_adata

# Compute the pearson correlation distance between the common features of the RNA and ATAC cells
mdata.uns["cross_rna+atac"] = cdist(merf_s_adata.X, atac_s_adata.X, metric="correlation")
mdata.uns["cross_keys"] = ["cross_rna+atac"]

mdata["rna"].layers["counts"] = mdata["rna"].X.copy()

# Log-normalize the scRNA counts
sc.pp.normalize_total(mdata["rna"], target_sum=10000.)
sc.pp.log1p(mdata["rna"])

# Perform PCA on the selected genes
sc.tl.pca(mdata["rna"], n_comps=100, zero_center=None)

import torch
torch.manual_seed(1792)
autoencoders = {"rna": scconfluence.unimodal.AutoEncoder(mdata["rna"],
                                                         modality="rna",
                                                         rep_in="X_pca",
                                                         rep_out="counts",
                                                         batch_key=None,
                                                         n_hidden=64,
                                                         n_latent=20,
                                                         type_loss="zinb"),
                "atac": scconfluence.unimodal.AutoEncoder(mdata["atac"],
                                                          modality="atac",
                                                          rep_in=None,
                                                          rep_out=None,
                                                          batch_key=None,
                                                          n_hidden=64,
                                                          n_latent=20,
                                                          type_loss="l2",
                                                          reconstruction_weight=5.)}

model = scconfluence.model.ScConfluence(mdata=mdata, unimodal_aes=autoencoders,
                                        mass=0.5, reach=0.3, iot_loss_weight=0.01, sinkhorn_loss_weight=0.1)
model.fit(save_path="Mop", use_cuda=False, max_epochs=1000)

mdata.obsm["latent"] = model.get_latent(use_cuda=False).loc[mdata.obs_names]

latent_data = mdata.obsm['latent']
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    latent_df = latent_data

latent_df.to_csv('./Reduced_feat/RF150/scCon.csv', index=False)