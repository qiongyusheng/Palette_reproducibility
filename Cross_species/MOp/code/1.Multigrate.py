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

dir = "./cross_species/"
n_batches = 8
rnas = []
atacs_h = []
atacs_m = []

rna_input = []
atac_h_input = []
atac_m_input = []

rna_t = [0,1,3,4,6,7]
atac_h_t = [0,2]
atac_m_t =[3,5]

path = ['human/SNARE_Lein','human/10XV3_Lein','human/Multiome_Zemke',
        'mouse/10XMultiome_Zemke','mouse/10XV3_Lein','mouse/10X',
       'Marmoset/SNARE_Lein','Marmoset/10XV3_Lein']

for batch in range(n_batches):
    
    if batch in rna_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/co.rds"))
        counts_rna = np.array(rds_data[None]).T
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch), counts_rna.shape[0])
        if batch in [0]:
            obs['modality'] = np.repeat('Multiome_human', counts_rna.shape[0])
        elif batch in [3]:
            obs['modality'] = np.repeat('Multiome_mouse', counts_rna.shape[0])
        else:
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
    
    if batch in atac_h_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/atac.rds"))
        counts_atac_h = np.array(rds_data[None]).T
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch), counts_atac_h.shape[0])
        if batch in [0]:
            obs['modality'] = np.repeat('Multiome_human', counts_atac_h.shape[0])
        else:
            obs['modality'] = np.repeat('atac_human', counts_atac_h.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_atac.shape[0])
        # np.repeat('atac', counts_atac.shape[0])
        atac_h_ad = ad.AnnData(counts_atac_h ,obs = obs)
        atac_h_ad.obs.index = rds_data[None].keys()
        sc.pp.normalize_total(atac_h_ad, target_sum=1e4)
        sc.pp.log1p(atac_h_ad)
        atac_h_ad.layers['log-norm'] = atac_h_ad.X.copy()
        atac_h_in = "log-norm"
        
    else:
        atac_h_ad = None
        atac_h_in = None
    

    if batch in atac_m_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/atac.rds"))
        counts_atac_m = np.array(rds_data[None]).T
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch), counts_atac_m.shape[0])
        if batch in [3]:
            obs['modality'] = np.repeat('Multiome_mouse', counts_atac_m.shape[0])
        else:
            obs['modality'] = np.repeat('atac_mouse', counts_atac_m.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_atac.shape[0])
        # np.repeat('atac', counts_atac.shape[0])
        atac_m_ad = ad.AnnData(counts_atac_m ,obs = obs)
        atac_m_ad.obs.index = rds_data[None].keys()
        sc.pp.normalize_total(atac_m_ad, target_sum=1e4)
        sc.pp.log1p(atac_m_ad)
        atac_m_ad.layers['log-norm'] = atac_m_ad.X.copy()
        atac_m_in = "log-norm"
        
    else:
        atac_m_ad = None
        atac_m_in = None
        
    rnas.append(rna_ad)
    rna_input.append(rna_in)
    
    atacs_h.append(atac_h_ad)
    atac_h_input.append(atac_h_in)
    
    atacs_m.append(atac_m_ad)
    atac_m_input.append(atac_m_in)

import multigrate as mtg

del rds_data
del counts_rna
del counts_atac_h
del counts_atac_m
import gc
gc.collect()

adata = mtg.data.organize_multiome_anndatas(
    adatas = [[rnas[0],rnas[1],rnas[2],rnas[3],
               rnas[4],rnas[5],rnas[6],rnas[7]], 
              
              [atacs_h[0],atacs_h[1],atacs_h[2],atacs_h[3],
               atacs_h[4],atacs_h[5],atacs_h[6],atacs_h[7]], 
              
              [atacs_m[0],atacs_m[1],atacs_m[2],atacs_m[3],
               atacs_m[4],atacs_m[5],atacs_m[6],atacs_m[7]]],
    
    layers = [[rna_input[0],rna_input[1],rna_input[2],rna_input[3],
               rna_input[4],rna_input[5],rna_input[6],rna_input[7]], 
              
              [atac_h_input[0],atac_h_input[1],atac_h_input[2],atac_h_input[3],
               atac_h_input[4],atac_h_input[5],atac_h_input[6],atac_h_input[7]],
              
              [atac_m_input[0],atac_m_input[1],atac_m_input[2],atac_m_input[3],
               atac_m_input[4],atac_m_input[5],atac_m_input[6],atac_m_input[7]]],
)

del rnas
del atacs_h
del atacs_m
gc.collect()

mtg.model.MultiVAE.setup_anndata(
    adata,
    categorical_covariate_keys=['modality', 'samplename'],
    rna_indices_end=400000, 
)

model = mtg.model.MultiVAE(
    adata,
    losses=['nb', 'mse', 'mse'],
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

latent_df.to_csv('./Cross_species/MOp/output/Multigrate.csv', index=False)
