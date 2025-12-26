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

dir = "./10x_Visium_human_tonsil/"
n_batches = 3
rnas = []
adts = []

rna_input = []
adt_input = []

rna_t = [0,1]
adt_t = [0,2]

path = ['slice1','slice2','slice3']

for batch in range(n_batches):
    
    if batch in rna_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/rna.rds"))
        counts_rna = np.array(rds_data[None]).T
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat(path[batch], counts_rna.shape[0])
        if batch in [0]:
            obs['modality'] = np.repeat('CITE', counts_rna.shape[0])
        else:
            obs['modality'] = np.repeat('RNA', counts_rna.shape[0])
        # np.repeat('rna', counts_rna.shape[0])
        rna_ad = ad.AnnData(counts_rna ,obs = obs)
        rna_ad.obs.index = rds_data[None].keys()
        rna_ad.layers["counts"] = rna_ad.X.copy()
        rna_in = "counts"
        
        # counts_rna = pyreadr.read_r(os.path.join(dir, 'rna' + str(batch + 1) + ".rds")).toarray().T
    else:
        rna_ad = None
        rna_in = None
    
    
    if batch in adt_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/adt.rds"))
        counts_protein = np.array(rds_data[None]).T
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat(path[batch], counts_protein.shape[0])
        if batch in [0]:
            obs['modality'] = np.repeat('CITE', counts_protein.shape[0])
        else:
            obs['modality'] = np.repeat('ADT', counts_protein.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_protein.shape[0])
        # np.repeat('adt', counts_protein.shape[0])
        adt_ad = ad.AnnData(counts_protein ,obs = obs)
        adt_ad.obs.index = rds_data[None].keys()
        muon.prot.pp.clr(adt_ad)
        adt_ad.layers['clr'] = adt_ad.X.copy()
        adt_in = "clr"
        
    else:
        adt_ad = None
        adt_in = None
        
    rnas.append(rna_ad)
    rna_input.append(rna_in)
    
    adts.append(adt_ad)
    adt_input.append(adt_in)

import multigrate as mtg

adata = mtg.data.organize_multiome_anndatas(
    adatas = [[rnas[0],rnas[1],rnas[2]], 
              [adts[0],adts[1],adts[2]]],
    layers = [[rna_input[0],rna_input[1],rna_input[2]], 
              [adt_input[0],adt_input[1],adt_input[2]]],
)

mtg.model.MultiVAE.setup_anndata(
    adata,
    categorical_covariate_keys=['modality', 'samplename'],
    rna_indices_end = 461
)

model = mtg.model.MultiVAE(
    adata,
    losses=['nb', 'mse'],
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

latent_df.to_csv('./10xVisium_tonsil/Multigrate.csv', index=False)