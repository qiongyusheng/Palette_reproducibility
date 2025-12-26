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
mdata = mu.MuData({"rna": rna_adata, "adt": adt_adata})

# Log-normalize scRNA counts
sc.pp.normalize_total(rna_s_adata, target_sum=10000.)
sc.pp.log1p(rna_s_adata)
# Use the Centered log ratio normalization for protein counts
mu.prot.pp.clr(adt_s_adata, axis=1)
sc.pp.scale(rna_s_adata)
sc.pp.scale(adt_s_adata)
mdata.uns["cross_rna+adt"] = cdist(rna_s_adata.X, adt_s_adata.X, metric="correlation")
mdata.uns["cross_keys"] = ["cross_rna+adt"]
mdata["rna"].layers["counts"] = mdata["rna"].X.copy()
# Log-normalize scRNA counts
sc.pp.normalize_total(mdata["rna"], target_sum=10000.)
sc.pp.log1p(mdata["rna"])

# Perform PCA on the selected genes
sc.tl.pca(mdata["rna"], n_comps=100, zero_center=None)

# Use the Centered log ratio normalization for protein counts
mu.prot.pp.clr(mdata["adt"], axis=1)
sc.tl.pca(mdata["adt"], n_comps=20, zero_center=None)

labels = pd.concat([
    pd.Series('RNA',  index=rna_adata.obs_names),
    pd.Series('ADT', index=adt_adata.obs_names),
])

# 对齐到 mdata 的观测索引再赋值
mdata.obs['modality'] = labels.reindex(mdata.obs_names).astype('category')
del labels

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
                "adt": scconfluence.unimodal.AutoEncoder(mdata["adt"],
                                                         modality="adt",
                                                         rep_in="X_pca",
                                                         rep_out=None,
                                                         batch_key=None,
                                                         n_hidden=64,
                                                         n_latent=20,
                                                         type_loss="l2",
                                                         reconstruction_weight=5.)}

model = scconfluence.model.ScConfluence(mdata=mdata, unimodal_aes=autoencoders,
                                        mass=0.5, reach=0.3, iot_loss_weight=0.01, sinkhorn_loss_weight=0.1)
model.fit(save_path="RNA_CyTOF", use_cuda=False, max_epochs=1000)
mdata.obsm["latent"] = model.get_latent(use_cuda=False).loc[mdata.obs_names]

latent_data = mdata.obsm['latent']
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    latent_df = latent_data

latent_df.to_csv('./scCon.csv', index=False)
