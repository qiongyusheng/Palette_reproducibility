import pyreadr
import sys, os
import numpy as np
from umap import UMAP
import time
import torch
import matplotlib.pyplot as plt
import pandas as pd  
import scipy.sparse as sp
from sklearn.decomposition import PCA
import scmomat.model as model
import scmomat.utils as utils
import scmomat.bmk as bmk
import scmomat.umap_batch as umap_batch
import scipy.io
from scipy import sparse
from scipy.stats import spearmanr, pearsonr
from scipy.spatial import distance_matrix
from scipy.sparse import coo_matrix
import anndata as ad
import scanpy as sc

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
def lsi(counts, n_components = 30):
    from sklearn.feature_extraction.text import TfidfTransformer
    from sklearn.decomposition import TruncatedSVD

    tfidf = TfidfTransformer(norm='l2', sublinear_tf=True)
    normed_count = tfidf.fit_transform(counts)

    # perform SVD on the sparse matrix
    lsi = TruncatedSVD(n_components=n_components + 1, random_state=42)
    lsi_r = lsi.fit_transform(normed_count)

    lsi.explained_variance_ratio_

    X_lsi = lsi_r[:, 1:]
    return X_lsi

dir = "./MERFISH_ATAC/"
n_batches = 2
counts_rnas = []
counts_atacs = []
labels = []
prec_labels = []
path = ['MERFISH','ATAC']

for batch in range(n_batches):
    if batch in [0,1]:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/co150_RF.rds"))
        counts_rna = np.array(rds_data[None]).T
        counts_rna = utils.preprocess(counts_rna, modality = "RNA", log = False)
        if batch in [1]:
            counts_rna = np.where(counts_rna > 0, 1, 0)
        # counts_rna = pyreadr.read_r(os.path.join(dir, 'rna' + str(batch + 1) + ".rds")).toarray().T
    else:
        counts_rna = None
    
    if batch in [1]:
        rds_data = scipy.io.mmread(os.path.join(dir, path[batch] + "/atac_midas.mtx")).transpose().tocsr()
        feat = pd.read_csv(os.path.join(dir, path[1] + "/feat.csv"))
        peak_name = np.array(feat.iloc[:, 0], dtype=str)
        adata = ad.AnnData(rds_data,dtype=np.float32)
        adata.var_names = peak_name
        
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=20000)
        adata_hvg = adata[:, adata.var.highly_variable].copy()
        regions_new = adata_hvg.var_names
        
        data = adata_hvg.X
        
        del adata
        del adata_hvg
        del feat
        del peak_name
        
        counts_atac = data.toarray()
        counts_atac = utils.preprocess(counts_atac, modality = "ATAC")
    else:
        counts_atac = None
        
    counts_rnas.append(counts_rna)
    counts_atacs.append(counts_atac)

counts = {"rna":counts_rnas, "atac": counts_atacs}

rds_data = pyreadr.read_r(os.path.join(dir, path[0] + "/co150_RF.rds"))
genes = rds_data[None].index
feats_name = {"rna": genes, "atac": regions_new}
counts["feats_name"] = feats_name

counts["nbatches"] = n_batches



model1 = model.scmomat_model(counts = counts,
                             device = device)
zs = []
for batch in range(n_batches):
    z = model1.softmax(model1.C_cells[str(batch)].cpu().detach()).numpy()
    z = z[:, :20]
    zs.append(z)

Z = np.concatenate(zs, axis=0)
latent_df = pd.DataFrame(Z)
latent_df.to_csv('./Reduced_feat/RF150/scMoMaT.csv',
                 index=False)
