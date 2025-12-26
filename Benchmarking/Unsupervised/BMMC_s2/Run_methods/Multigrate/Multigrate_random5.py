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
import scipy.io as sio

import warnings
warnings.filterwarnings("ignore")


dir = "./BMMC_CITE_Multiome/"
n_batches = 10

rnas = []
atacs = []
adts = []

rna_input = []
atac_input = []
adt_input = []

rna_t = [0,None,
         1,
         2,None,
         3,None,
         4,None,
         5]

atac_t = [None,None,
          None,None,
          None,
          None,3,
          None,4,
          5]

adt_t =[None,0,
        1,
        None,2,
        None,None,
        None,
        None,None]

path = ['CITE/s1d1_cite','CITE/s1d2_cite','CITE/s1d3_cite','Multiome/s1d1_multi','Multiome/s1d2_multi','Multiome/s1d3_multi']


s_name = ['batch1','batch2','batch3','batch4','batch5',
          'batch6','batch7','batch8','batch9','batch10']


modal_name = ['rna','adt',
              'cite',
              'rna','adt',
              'rna','atac',
              'rna','atac',
              'multiome']


for batch in range(n_batches):
    
    if rna_t[batch] in [0,1,2,3,4,5]:
        rds_data = pyreadr.read_r(os.path.join(dir, path[rna_t[batch]] + "/rna.rds"))
        # counts_rna = np.array(rds_data[None]).T
        counts_rna = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat(s_name[batch], counts_rna.shape[0])
        obs['modality'] = np.repeat(modal_name[batch], counts_rna.shape[0])
        rna_ad = ad.AnnData(counts_rna ,obs = obs)
        rna_ad.obs.index = rds_data[None].keys()
        if rna_t[batch] in [0,2,4,5]:
            rna_ad.obs.index = [f"{idx}_RNA" for idx in rna_ad.obs.index]
        rna_ad.layers["counts"] = rna_ad.X.copy()
        rna_in = "counts"
        
        # counts_rna = pyreadr.read_r(os.path.join(dir, 'rna' + str(batch + 1) + ".rds")).toarray().T
    else:
        rna_ad = None
        rna_in = None
    
    if atac_t[batch] in [3,4,5]:
        rds_data = pyreadr.read_r(os.path.join(dir, path[atac_t[batch]] + "/atac.rds"))
        # counts_atac = np.array(rds_data[None]).T
        counts_atac = csr_matrix(np.array(rds_data[None]).T)
        
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat(s_name[batch], counts_atac.shape[0])
        obs['modality'] = np.repeat(modal_name[batch], counts_atac.shape[0])    
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_atac.shape[0])
        # np.repeat('atac', counts_atac.shape[0])
        atac_ad = ad.AnnData(counts_atac ,obs = obs)
        atac_ad.obs.index = rds_data[None].keys()
        if atac_t[batch] in [0,2,4,5]:
            atac_ad.obs.index = [f"{idx}_ATAC" for idx in atac_ad.obs.index]
        sc.pp.normalize_total(atac_ad, target_sum=1e4)
        sc.pp.log1p(atac_ad)
        atac_ad.layers['log-norm'] = atac_ad.X.copy()
        atac_in = "log-norm"
        
    else:
        atac_ad = None
        atac_in = None
    
    if adt_t[batch] in [0,1,2]:
        rds_data = pyreadr.read_r(os.path.join(dir, path[adt_t[batch]] + "/adt.rds"))
        # counts_protein = np.array(rds_data[None]).T
        counts_protein = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat(s_name[batch], counts_protein.shape[0])
        obs['modality'] = np.repeat(modal_name[batch], counts_protein.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_protein.shape[0])
        # np.repeat('adt', counts_protein.shape[0])
        adt_ad = ad.AnnData(counts_protein ,obs = obs)
        adt_ad.obs.index = rds_data[None].keys()
        if adt_t[batch] in [0,2,4,5]:
            adt_ad.obs.index = [f"{idx}_ADT" for idx in adt_ad.obs.index]
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

adata_combined = sc.concat([atacs[6], atacs[8], atacs[9]], join='inner')
sc.pp.highly_variable_genes(adata_combined, n_top_genes=45000, batch_key='samplename')
adata_hvg = adata_combined[:, adata_combined.var.highly_variable].copy()
adata_hvg6 = adata_hvg[adata_combined.obs['samplename'] == s_name[6],:]
adata_hvg8 = adata_hvg[adata_combined.obs['samplename'] == s_name[8],:]
adata_hvg9 = adata_hvg[adata_combined.obs['samplename'] == s_name[9],:]

atacs[6] = adata_hvg6
atacs[8] = adata_hvg8
atacs[9] = adata_hvg9

import multigrate as mtg
del rds_data
del counts_protein
del counts_atac
del counts_rna
import gc
gc.collect()

adata = mtg.data.organize_multiome_anndatas(
    adatas = [[rnas[0],rnas[1],rnas[2],rnas[3],rnas[4],
               rnas[5],rnas[6],rnas[7],rnas[8],rnas[9]], 
              
              [atacs[0],atacs[1],atacs[2],atacs[3],atacs[4],
               atacs[5],atacs[6],atacs[7],atacs[8],atacs[9]], 
              
              [adts[0],adts[1],adts[2],adts[3],adts[4],
               adts[5],adts[6],adts[7],adts[8],adts[9]]],
    
    layers = [[rna_input[0],rna_input[1],rna_input[2],rna_input[3],rna_input[4],
               rna_input[5],rna_input[6],rna_input[7],rna_input[8],rna_input[9]], 
              
              [atac_input[0],atac_input[1],atac_input[2],atac_input[3],atac_input[4],
               atac_input[5],atac_input[6],atac_input[7],atac_input[8],atac_input[9]], 
              
              [adt_input[0],adt_input[1],adt_input[2],adt_input[3],adt_input[4],
               adt_input[5],adt_input[6],adt_input[7],adt_input[8],adt_input[9]]],
)

del rnas
del atacs
del adts
gc.collect()

mtg.model.MultiVAE.setup_anndata(
    adata,
    categorical_covariate_keys=['modality', 'samplename'],
    rna_indices_end = 1000000
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

# 保存为 CSV 文件
latent_df.to_csv('./Benchmarking/Unsupervised/BMMC_s2/output/Multigrate_26.csv', index=False)