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

dir = "./cross_species/"
n_batches = 8
counts_rnas = []
counts_atacs_h = []
counts_atacs_m = []
labels = []
prec_labels = []
path = ['human/SNARE_Lein','human/10XV3_Lein','human/Multiome_Zemke',
        'mouse/10XMultiome_Zemke','mouse/10XV3_Lein','mouse/10X',
       'Marmoset/SNARE_Lein','Marmoset/10XV3_Lein']

for batch in range(n_batches):
    if batch in [0,1,3,4,6,7]:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/co.rds"))
        counts_rna = np.array(rds_data[None]).T
        counts_rna = utils.preprocess(counts_rna, modality = "RNA", log = False)
        # counts_rna = pyreadr.read_r(os.path.join(dir, 'rna' + str(batch + 1) + ".rds")).toarray().T
    else:
        counts_rna = None
    
    if batch in [0,2]:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/atac.rds"))
        counts_atac_h = np.array(rds_data[None]).T
        counts_atac_h = utils.preprocess(counts_atac_h, modality = "ATAC")
    else:
        counts_atac_h = None
        
    
    if batch in [3,5]:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/atac.rds"))
        counts_atac_m = np.array(rds_data[None]).T
        counts_atac_m = utils.preprocess(counts_atac_m, modality = "ATAC")
    else:
        counts_atac_m = None
    

        
    counts_rnas.append(counts_rna)
    counts_atacs_h.append(counts_atac_h)
    counts_atacs_m.append(counts_atac_m)

counts = {"rna":counts_rnas, "atac_h": counts_atacs_h, "atac_m": counts_atacs_m}

for modality, data_list in counts.items():
    print(f"Modality: {modality}")
    for idx, data in enumerate(data_list):
        if data is not None:
            print(f"Data {idx} shape: {data.shape}")
        else:
            print(f"Data {idx} is None")
    print()

rds_data = pyreadr.read_r(os.path.join(dir, path[0] + "/co.rds"))
genes = rds_data[None].index
rds_data = pyreadr.read_r(os.path.join(dir, path[0] + "/atac.rds"))
regions_h = rds_data[None].index
rds_data = pyreadr.read_r(os.path.join(dir, path[5] + "/atac.rds"))
regions_m = rds_data[None].index
feats_name = {"rna": genes, "atac_h": regions_h, "atac_m": regions_m}
counts["feats_name"] = feats_name

counts["nbatches"] = n_batches

del rds_data
del counts_atac_h
del counts_atac_m
del counts_rna
del counts_rnas
del counts_atacs_h
del counts_atacs_m


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
latent_df.to_csv('./Cross_species/MOp/output/scMoMaT.csv', index=False)
