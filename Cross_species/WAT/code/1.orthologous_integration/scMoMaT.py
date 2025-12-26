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


dir = "/home/server/sqy/data/mosaic/input/human_mouse_WAT/"
n_batches = 36
counts_rnas = []

labels = []
prec_labels = []

path = ['human/Hs012','human/Hs013','human/Hs266','human/Hs001','human/Hs002',
        'human/Hs004','human/Hs253','human/Hs254','human/Hs256','human/Hs255',
        'human/Hs009','human/Hs010','human/Hs011','human/Hs235','human/Hs236',
        'human/Hs237','human/Hs238','human/Hs239','human/Hs240','human/Hs242',
        'human/Hs248','human/Hs249',
        'mouse/HFD06.F','mouse/NCD01.F','mouse/HFD07.F','mouse/NCD02.F',
        'mouse/HFD01','mouse/NCD06','mouse/HFD02','mouse/HFD03','mouse/HFD04',
        'mouse/HFD05','mouse/NCD09','mouse/NCD10','mouse/NCD07','mouse/NCD08']

for batch in range(n_batches):
    rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/homo.rds"))
    counts_rna = np.array(rds_data[None]).T
    counts_rna = utils.preprocess(counts_rna, modality = "RNA", log = False)
        
    counts_rnas.append(counts_rna)
counts = {"rna":counts_rnas}

del counts_rnas
del rds_data

for modality, data_list in counts.items():
    print(f"Modality: {modality}")
    for idx, data in enumerate(data_list):
        if data is not None:
            print(f"Data {idx} shape: {data.shape}")
        else:
            print(f"Data {idx} is None")
    print()

rds_data = pyreadr.read_r(os.path.join(dir, path[0] + "/homo.rds"))
genes = rds_data[None].index

feats_name = {"rna": genes}
counts["feats_name"] = feats_name

counts["nbatches"] = n_batches

del rds_data

model1 = model.scmomat_model(counts = counts, 
                             device = device)
losses1 = model1.train_func()

zs = []
for batch in range(n_batches):
    z = model1.softmax(model1.C_cells[str(batch)].cpu().detach()).numpy()
    z = z[:, :20]
    zs.append(z)

Z = np.concatenate(zs, axis=0)
latent_df = pd.DataFrame(Z)
latent_df.to_csv('./scMoMaT_homo.csv', index=False)

