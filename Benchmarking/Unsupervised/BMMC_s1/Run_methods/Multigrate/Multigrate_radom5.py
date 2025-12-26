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
from scipy.sparse import csr_matrix
from scipy import io
import scipy.sparse as sp
import warnings

warnings.filterwarnings("ignore")
dir = "./BMMC_CITE_Multiome/"
n_batches = 5
rnas = []
atacs = []
adts = []

rna_input = []
atac_input = []
adt_input = []

rna_t = [0,1,2,3,4]
atac_t = [2,3,4]
adt_t =[0,1]

path = ['CITE/s1d1_cite','CITE/s1d3_cite',
        'Multiome/s1d1_multi','Multiome/s1d2_multi','Multiome/s1d3_multi']

for batch in range(n_batches):
    
    if batch in rna_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/rna.rds"))
        # counts_rna = np.array(rds_data[None]).T
        counts_rna = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch), counts_rna.shape[0])
        if batch in [0,1,2]:
            obs['modality'] = np.repeat('CITE', counts_rna.shape[0])
        else:
            obs['modality'] = np.repeat('Multiome', counts_rna.shape[0])
        # np.repeat('rna', counts_rna.shape[0])
        rna_ad = ad.AnnData(counts_rna ,obs = obs)
        rna_ad.obs.index = rds_data[None].keys()
        rna_ad.layers["counts"] = rna_ad.X.copy()
        rna_in = "counts"
        
        # counts_rna = pyreadr.read_r(os.path.join(dir, 'rna' + str(batch + 1) + ".rds")).toarray().T
    else:
        rna_ad = None
        rna_in = None
    
    if batch in atac_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/atac.rds"))
        # counts_atac = np.array(rds_data[None]).T
        counts_atac = csr_matrix(np.array(rds_data[None]).T)
        
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch), counts_atac.shape[0])
        obs['modality'] = np.repeat('Multiome', counts_atac.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_atac.shape[0])
        # np.repeat('atac', counts_atac.shape[0])
        atac_ad = ad.AnnData(counts_atac ,obs = obs)
        atac_ad.obs.index = rds_data[None].keys()
        sc.pp.normalize_total(atac_ad, target_sum=1e4)
        sc.pp.log1p(atac_ad)
        atac_ad.layers['log-norm'] = atac_ad.X.copy()
        atac_in = "log-norm"
        
    else:
        atac_ad = None
        atac_in = None
    
    if batch in adt_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/adt.rds"))
        # counts_protein = np.array(rds_data[None]).T
        counts_protein = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch), counts_protein.shape[0])
        obs['modality'] = np.repeat('CITE', counts_protein.shape[0])
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
    
    atacs.append(atac_ad)
    atac_input.append(atac_in)
    
    adts.append(adt_ad)
    adt_input.append(adt_in)

import multigrate as mtg
del rds_data
del counts_protein
del counts_atac
del counts_rna
import gc
gc.collect()

adata = mtg.data.organize_multiome_anndatas(
    adatas = [[rnas[0],rnas[1],rnas[2],rnas[3],rnas[4]], 
              [atacs[0],atacs[1],atacs[2],atacs[3],atacs[4]], 
              [adts[0],adts[1],adts[2],adts[3],adts[4]]],
    layers = [[rna_input[0],rna_input[1],rna_input[2],rna_input[3],rna_input[4]], 
              [atac_input[0],atac_input[1],atac_input[2],atac_input[3],atac_input[4]], 
              [adt_input[0],adt_input[1],adt_input[2],adt_input[3],adt_input[4]]],
)

del rnas
del atacs
del adts
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

latent_df.to_csv('./Benchmarking/Unsupervised/BMMC_s1/output/Multigrate_drop2.csv', index=False)