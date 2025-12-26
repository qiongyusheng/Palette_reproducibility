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
import scipy.io as sio
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

dir = "./Human_Retina/"
n_batches = 8
counts_rnas = []
counts_atacs = []
counts_proteins = []
labels = []
prec_labels = []
path = ['LGS1OD',
        'LGS1OS',
        'LGS2OD',
        'LGS2OS',
        'LGS3OD',
        'LGS3OS',
        'LVG1OD',
        'LVG1OS']

for batch in range(n_batches):
    if batch in [0,2,4,5,6]:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/rna.rds"))
        counts_rna = np.array(rds_data[None]).T
        counts_rna = utils.preprocess(counts_rna, modality = "RNA", log = False)
        if batch in [6]:
            genes = rds_data[None].index
    else:
        counts_rna = None
    
    if batch in [1,3,5,6,7]:
        matrix = sio.mmread(os.path.join(dir, path[batch] + "/atac.mtx"))
        counts_atac = matrix.transpose().tocsr()
        if batch in [6]:
            peak_name = pd.read_csv(os.path.join(dir, path[batch] + "/peak.csv"))
            regions = pd.Index(peak_name.iloc[:, 0], dtype='object', name='rownames')
    else:
        counts_atac = None
        
    counts_rnas.append(counts_rna)
    counts_atacs.append(counts_atac)

import scipy.sparse as sp
tmp = sp.vstack([counts_atacs[1], counts_atacs[3], counts_atacs[5], counts_atacs[6], counts_atacs[7]])
import anndata as ad
import scanpy as sc
adata = ad.AnnData(tmp,dtype=np.float32)
adata.obs['batch'] = ["batch1"] * counts_atacs[1].shape[0] + ["batch2"] * counts_atacs[3].shape[0] + ["batch3"] * counts_atacs[5].shape[0] + ["batch4"] * counts_atacs[6].shape[0] + ["batch5"] * counts_atacs[7].shape[0]
# adata.layers['counts'] = adata.X.copy()
adata.var_names = peak_name
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=50000, batch_key='batch')

adata_hvg = adata[:, adata.var.highly_variable].copy()
regions_new = adata_hvg.var_names

adata_hvg0 = adata_hvg[adata.obs['batch'] == 'batch1',:]
hvg0_d = np.array(adata_hvg0.X.todense())
counts_atac0 = utils.preprocess(hvg0_d, modality = "ATAC")
counts_atacs[1] = counts_atac0
del adata_hvg0
del counts_atac0
del hvg0_d

adata_hvg1 = adata_hvg[adata.obs['batch'] == 'batch2',:]
hvg1_d = np.array(adata_hvg1.X.todense())
counts_atac1 = utils.preprocess(hvg1_d, modality = "ATAC")
counts_atacs[3] = counts_atac1
del adata_hvg1
del counts_atac1
del hvg1_d

adata_hvg3 = adata_hvg[adata.obs['batch'] == 'batch3',:]
hvg3_d = np.array(adata_hvg3.X.todense())
counts_atac3 = utils.preprocess(hvg3_d, modality = "ATAC")
counts_atacs[5] = counts_atac3
del adata_hvg3
del counts_atac3
del hvg3_d

adata_hvg6 = adata_hvg[adata.obs['batch'] == 'batch4',:]
hvg6_d = np.array(adata_hvg6.X.todense())
counts_atac6 = utils.preprocess(hvg6_d, modality = "ATAC")
counts_atacs[6] = counts_atac6
del adata_hvg6
del counts_atac6
del hvg6_d

adata_hvg5 = adata_hvg[adata.obs['batch'] == 'batch5',:]
hvg5_d = np.array(adata_hvg5.X.todense())
counts_atac5 = utils.preprocess(hvg5_d, modality = "ATAC")
counts_atacs[7] = counts_atac5
del adata_hvg5
del counts_atac5
del hvg5_d

del tmp
del adata
del adata_hvg

counts = {"rna":counts_rnas, "atac": counts_atacs}
feats_name = {"rna": genes, "atac": regions_new}
counts["feats_name"] = feats_name
counts["nbatches"] = n_batches

del counts_rnas
del counts_atacs
del counts_rna
del counts_atac
del rds_data
import gc
gc.collect()

model1 = model.scmomat_model(counts = counts, 
                             device = device)
zs = []
for batch in range(n_batches):
    z = model1.softmax(model1.C_cells[str(batch)].cpu().detach()).numpy()
    z = z[:, :20]
    zs.append(z)

Z = np.concatenate(zs, axis=0)
latent_df = pd.DataFrame(Z)
latent_df.to_csv('./Benchmarking/Unsupervised/Retina/output/scMoMaT_s4.csv', index=False)