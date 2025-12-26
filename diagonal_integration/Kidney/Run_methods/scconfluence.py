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

atac = sc.read('./Muto-2021-ATAC.h5ad')
rna = sc.read('./Muto-2021-RNA.h5ad')

atac.obs_names = atac.obs_names + '_ATAC'
rna.obs_names = rna.obs_names + '_RNA'

rds_rna = pyreadr.read_r('./RNA/X.rds')
feat = rds_rna[None].index
feat = list(feat.astype(str)) if hasattr(feat, "astype") else list(map(str, feat))
gene_ids = pd.Index(rna.var['gene_ids'].astype(str))   # 关键：转成 Index
idx = gene_ids.get_indexer(feat)   # Index 方法
present_mask = idx >= 0
idx_keep = idx[present_mask]
rna = rna[:, idx_keep].copy()
del rds_rna
del feat
del gene_ids
del present_mask
del idx
del idx_keep
mdata = mu.MuData({"rna": rna, "atac": atac})
labels = pd.concat([
    pd.Series('RNA',  index=rna.obs_names),
    pd.Series('ATAC', index=atac.obs_names),
])

mdata.obs['modality'] = labels.reindex(mdata.obs_names).astype('category')
del labels

atac_cm = pyreadr.read_r('./ATAC/co.rds')
rna_cm = pyreadr.read_r('./RNA/co.rds')
cm_features_atac = ad.AnnData(csr_matrix(np.array(atac_cm[None]).T))
cm_features_atac.obs_names = atac_cm[None].columns
cm_features_rna = ad.AnnData(csr_matrix(np.array(rna_cm[None]).T))
cm_features_rna.obs_names = rna_cm[None].columns + '_RNA'
sc.pp.normalize_total(cm_features_rna, target_sum=10000.)
sc.pp.log1p(cm_features_rna)

# Log-normalize scATAC gene activities
sc.pp.normalize_total(cm_features_atac, target_sum=10000.)
sc.pp.log1p(cm_features_atac)

# Select highly variable genes for both modalities using the more reliable scRNA counts
cm_hvg_genes = sc.pp.highly_variable_genes(cm_features_rna, n_top_genes=1565, subset=False, inplace=False)
cm_features_rna = cm_features_rna[:, cm_hvg_genes["highly_variable"]].copy()
cm_features_atac = cm_features_atac[:, cm_hvg_genes["highly_variable"]].copy()

sc.pp.scale(cm_features_rna)
sc.pp.scale(cm_features_atac)

# Compute the pearson correlation distance between the common features of the RNA and ATAC cells
mdata.uns["cross_rna+atac"] = cdist(cm_features_rna.X, cm_features_atac.X, metric="correlation")
mdata.uns["cross_keys"] = ["cross_rna+atac"]

mdata["rna"].layers["counts"] = mdata["rna"].X.copy()

# Log-normalize the scRNA counts
sc.pp.normalize_total(mdata["rna"], target_sum=10000.)
sc.pp.log1p(mdata["rna"])

# Perform PCA on the selected genes
sc.tl.pca(mdata["rna"], n_comps=100, zero_center=None)

mu.atac.pp.tfidf(mdata["atac"], log_tf=True, log_idf=True)
sc.tl.pca(mdata["atac"], n_comps=100, zero_center=None)

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
                                                          rep_in="X_pca",
                                                          rep_out=None,
                                                          batch_key=None,
                                                          n_hidden=64,
                                                          n_latent=20,
                                                          type_loss="l2",
                                                          reconstruction_weight=5.)}

model = scconfluence.model.ScConfluence(mdata=mdata, unimodal_aes=autoencoders,
                                        mass=0.5, reach=0.3, iot_loss_weight=0.01, sinkhorn_loss_weight=0.1)
model.fit(save_path="kidney_rna_atac", use_cuda=False, max_epochs=1000)

mdata.obsm["latent"] = model.get_latent(use_cuda=False).loc[mdata.obs_names]

latent_data = mdata.obsm['latent']
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    latent_df = latent_data

latent_df.to_csv('./kidney/scCon.csv', index=False)
