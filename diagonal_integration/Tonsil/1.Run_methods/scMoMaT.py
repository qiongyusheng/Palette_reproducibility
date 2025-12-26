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

dir = "./"
n_batches = 2
counts_rnas = []
counts_proteins = []
labels = []
prec_labels = []
path = ['RNA','CODEX']

for batch in range(n_batches):
    if batch in [0]:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/X.rds"))
        counts_rna = np.array(rds_data[None]).T
        counts_rna = utils.preprocess(counts_rna, modality = "RNA", log = False)
        # counts_rna = pyreadr.read_r(os.path.join(dir, 'rna' + str(batch + 1) + ".rds")).toarray().T
    else:
        counts_rna = None

    if batch in [0,1]:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/co_all.rds"))
        counts_protein = np.array(rds_data[None]).T
        # counts_protein = utils.preprocess(counts_protein, modality = "RNA", log = True)
    else:
        counts_protein = None 
        
    counts_rnas.append(counts_rna)
    counts_proteins.append(counts_protein)

counts = {"rna":counts_rnas,  "protein": counts_proteins}

rds_data = pyreadr.read_r(os.path.join(dir, path[0] + "/rna.rds"))
genes = rds_data[None].index
rds_data = pyreadr.read_r(os.path.join(dir, path[0] + "/co_all.rds"))
proteins = rds_data[None].index
feats_name = {"rna": genes,  "protein": proteins}
counts["feats_name"] = feats_name
counts["nbatches"] = n_batches

model1 = model.scmomat_model(counts = counts, 
                             device = device)
zs = []
for batch in range(n_batches):
    z = model1.softmax(model1.C_cells[str(batch)].cpu().detach()).numpy()
    z = z
    zs.append(z)

Z = np.concatenate(zs, axis=0)
latent_df = pd.DataFrame(Z)

latent_df.to_csv('./scMoMaT.csv', index=False)