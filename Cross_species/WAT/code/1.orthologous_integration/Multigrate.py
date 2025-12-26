import scarches as sca
import scanpy as sc
import anndata as ad
import numpy as np
import muon
import gdown
import json
import os
import pyreadr
import pandas as pd

import warnings
warnings.filterwarnings("ignore")

dir = "./human_mouse_WAT/"
n_batches = 36
rnas = []

rna_input = []

rna_t = range(n_batches)

path = ['human/Hs012','human/Hs013','human/Hs266','human/Hs001','human/Hs002',
        'human/Hs004','human/Hs253','human/Hs254','human/Hs256','human/Hs255',
        'human/Hs009','human/Hs010','human/Hs011','human/Hs235','human/Hs236',
        'human/Hs237','human/Hs238','human/Hs239','human/Hs240','human/Hs242',
        'human/Hs248','human/Hs249',
        'mouse/HFD06.F','mouse/NCD01.F','mouse/HFD07.F','mouse/NCD02.F',
        'mouse/HFD01','mouse/NCD06','mouse/HFD02','mouse/HFD03','mouse/HFD04',
        'mouse/HFD05','mouse/NCD09','mouse/NCD10','mouse/NCD07','mouse/NCD08']

for batch in range(n_batches):
    
    if batch in rna_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/homo.rds"))
        counts_rna = np.array(rds_data[None]).T
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat(path[batch], counts_rna.shape[0])
        obs['modality'] = np.repeat('rna', counts_rna.shape[0])
        # np.repeat('rna', counts_rna.shape[0])
        rna_ad = ad.AnnData(counts_rna ,obs = obs)
        rna_ad.obs.index = rds_data[None].keys()
        rna_ad.layers["counts"] = rna_ad.X.copy()
        rna_in = "counts"
        
        # counts_rna = pyreadr.read_r(os.path.join(dir, 'rna' + str(batch + 1) + ".rds")).toarray().T
    else:
        rna_ad = None
        rna_in = None

    rnas.append(rna_ad)
    rna_input.append(rna_in)

import multigrate as mtg

adata = mtg.data.organize_multiome_anndatas(
    adatas = [[rnas[0],rnas[1],rnas[2],
              rnas[3],rnas[4],rnas[5],
              rnas[6],rnas[7],rnas[8],
              rnas[9],rnas[10],rnas[11],
              rnas[12],rnas[13],rnas[14],
              rnas[15],rnas[16],rnas[17],
              rnas[18],rnas[19],rnas[20],
              rnas[21],rnas[22],rnas[23],
              rnas[24],rnas[25],rnas[26],
              rnas[27],rnas[28],rnas[29],
              rnas[30],rnas[31],rnas[32],
              rnas[33],rnas[34],rnas[35]],],
    layers = [[rna_input[0],rna_input[1],rna_input[2],
              rna_input[3],rna_input[4],rna_input[5],
              rna_input[6],rna_input[7],rna_input[8],
              rna_input[9],rna_input[10],rna_input[11],
              rna_input[12],rna_input[13],rna_input[14],
              rna_input[15],rna_input[16],rna_input[17],
              rna_input[18],rna_input[19],rna_input[20],
              rna_input[21],rna_input[22],rna_input[23],
              rna_input[24],rna_input[25],rna_input[26],
              rna_input[27],rna_input[28],rna_input[29],
              rna_input[30],rna_input[31],rna_input[32],
              rna_input[33],rna_input[34],rna_input[35]]],
)

del rnas
del rds_data
del counts_rna
import gc
gc.collect()
mtg.model.MultiVAE.setup_anndata(
    adata,
    categorical_covariate_keys=['modality', 'samplename'],
    rna_indices_end = 5000
)

model = mtg.model.MultiVAE(
    adata,
    losses=['nb'],
    loss_coefs={'kl': 1e-3,
               'integ': 3000,
               },
    integrate_on='modality',
    mmd='marginal',z_dim = 20
)

model.train()
model.get_latent_representation()

latent_data = adata.obsm['latent']
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    latent_df = latent_data

latent_df.to_csv('./Multigrate_homo.csv', index=False)
