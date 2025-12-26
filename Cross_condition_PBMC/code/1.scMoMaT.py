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
path = ['PBMC_ASAP/control_asap','PBMC_ASAP/stim_asap',
          'PBMC_CITE/control_cite','PBMC_CITE/stim_cite']

n_batches = 4
counts_rnas = []
counts_atacs = []
counts_proteins = []
labels = []
prec_labels = []

for batch in range(n_batches):
    rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/adt.rds"))
    counts_protein = np.array(rds_data[None]).T
    counts_protein = utils.preprocess(counts_protein, modality = "RNA", log = True)
    
    if batch in [2,3]:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/rna.rds"))
        counts_rna = np.array(rds_data[None]).T
        counts_rna = utils.preprocess(counts_rna, modality = "RNA", log = False)
        # counts_rna = pyreadr.read_r(os.path.join(dir, 'rna' + str(batch + 1) + ".rds")).toarray().T
    else:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/pseudo_rna.rds"))
        counts_rna = np.array(rds_data[None]).T
    
    if batch in [0,1]:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/atac.rds"))
        counts_atac = np.array(rds_data[None]).T
        counts_atac = utils.preprocess(counts_atac, modality = "ATAC")
    else:
        counts_atac = None
    
        
    counts_rnas.append(counts_rna)
    counts_atacs.append(counts_atac)
    counts_proteins.append(counts_protein)

counts = {"rna":counts_rnas, "atac": counts_atacs, "protein": counts_proteins}

for modality, data_list in counts.items():
    print(f"Modality: {modality}")
    for idx, data in enumerate(data_list):
        if data is not None:
            print(f"Data {idx} shape: {data.shape}")
        else:
            print(f"Data {idx} is None")
    print()

rds_data = pyreadr.read_r(os.path.join(dir, path[3] + "/rna.rds"))
genes = rds_data[None].index
rds_data = pyreadr.read_r(os.path.join(dir, path[0] + "/adt.rds"))
proteins = rds_data[None].index
rds_data = pyreadr.read_r(os.path.join(dir, path[0] + "/atac.rds"))
regions = rds_data[None].index
feats_name = {"rna": genes, "atac": regions, "protein": proteins}
counts["feats_name"] = feats_name

counts["nbatches"] = n_batches

model1 = model.scmomat_model(counts = counts, 
                             device = device)
losses1 = model1.train_func()


zs = []
for batch in range(n_batches):
    z = model1.softmax(model1.C_cells[str(batch)].cpu().detach()).numpy()
    z = z[:, :20]
    zs.append(z)
    
Z = np.concatenate(zs, axis = 0)

latent_df = pd.DataFrame(Z)
latent_df.to_csv('./Cross_condition_PBMC/output/Embedding/scMoMaT.csv', index=False)
