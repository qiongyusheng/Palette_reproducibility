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


dir = "./Human_Retina/"
n_batches = 8
rnas = []
atacs = []

rna_input = []
atac_input = []

rna_t = [2,3,4,5,7]
atac_t = [0,1,3,5,6]

path = ['LGS1OD',
        'LGS1OS',
        'LGS2OD',
        'LGS2OS',
        'LGS3OD',
        'LGS3OS',
        'LVG1OD',
        'LVG1OS']

for batch in range(n_batches):
    
    if batch in rna_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/rna.rds"))
        # counts_rna = np.array(rds_data[None]).T
        counts_rna = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch), counts_rna.shape[0])
        if batch in [3,5]:
            obs['modality'] = np.repeat('Multiome', counts_rna.shape[0])
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
    
    if batch in atac_t:
        # rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/atac.rds"))
        # counts_atac = np.array(rds_data[None]).T
        matrix = sio.mmread(os.path.join(dir, path[batch] + "/atac.mtx"))
        counts_atac = matrix.transpose().tocsr()
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch), counts_atac.shape[0])
        # obs['modality'] = np.repeat('Multiome', counts_atac.shape[0])
        if batch in [3,5]:
            obs['modality'] = np.repeat('Multiome', counts_atac.shape[0])
        else:
            obs['modality'] = np.repeat('ATAC', counts_atac.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_atac.shape[0])
        # np.repeat('atac', counts_atac.shape[0])
        atac_ad = ad.AnnData(counts_atac ,obs = obs)
        barcode = pd.read_csv(os.path.join(dir, path[batch] + "/bcd.csv"))
        atac_ad.obs.index = pd.Index(barcode.iloc[:, 0], dtype='object')
        # atac_ad.obs.index = rds_data[None].keys()
        sc.pp.normalize_total(atac_ad, target_sum=1e4)
        sc.pp.log1p(atac_ad)
        atac_ad.layers['log-norm'] = atac_ad.X.copy()
        atac_in = "log-norm"
        
    else:
        atac_ad = None
        atac_in = None
    
        
    rnas.append(rna_ad)
    rna_input.append(rna_in)
    
    atacs.append(atac_ad)
    atac_input.append(atac_in)

adata_combined = ad.concat([atacs[0],atacs[1],atacs[3],atacs[5],atacs[6]], axis=0)
peak = pd.read_csv(os.path.join(dir, path[0] + "/peak.csv"))
adata_combined.var_names = pd.Index(peak.iloc[:, 0], dtype='object')
sc.pp.highly_variable_genes(adata_combined, n_top_genes=50000, batch_key='samplename')
adata_hvg = adata_combined[:, adata_combined.var.highly_variable].copy()
peak = pd.DataFrame(adata_hvg.var_names)

adata_hvg0 = adata_hvg[adata_combined.obs['samplename'] == 'batch0',:]
adata_hvg1 = adata_hvg[adata_combined.obs['samplename'] == 'batch1',:]
adata_hvg5 = adata_hvg[adata_combined.obs['samplename'] == 'batch3',:]
adata_hvg6 = adata_hvg[adata_combined.obs['samplename'] == 'batch5',:]
adata_hvg7 = adata_hvg[adata_combined.obs['samplename'] == 'batch6',:]

atacs[0] = adata_hvg0
atacs[1] = adata_hvg1
atacs[3] = adata_hvg5
atacs[5] = adata_hvg6
atacs[6] = adata_hvg7

del adata_combined
del adata_hvg
del adata_hvg0
del adata_hvg1
del adata_hvg5
del adata_hvg6
del adata_hvg7

import multigrate as mtg
del rds_data
del matrix
del counts_atac
del counts_rna
import gc
gc.collect()

adata = mtg.data.organize_multiome_anndatas(
    adatas = [[rnas[0],rnas[1],rnas[2],rnas[3],rnas[4],
               rnas[5],rnas[6],rnas[7]], 
              
              [atacs[0],atacs[1],atacs[2],atacs[3],atacs[4],
               atacs[5],atacs[6],atacs[7]]],
    
    layers = [[rna_input[0],rna_input[1],rna_input[2],rna_input[3],rna_input[4],
               rna_input[5],rna_input[6],rna_input[7]], 
              
              [atac_input[0],atac_input[1],atac_input[2],atac_input[3],atac_input[4],
               atac_input[5],atac_input[6],atac_input[7]]],
)
del rnas
del atacs
gc.collect()

mtg.model.MultiVAE.setup_anndata(
    adata,
    categorical_covariate_keys=['modality', 'samplename'],
    rna_indices_end=400000, 
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

# 保存为 CSV 文件
latent_df.to_csv('./Benchmarking/Unsupervised/Retina/output/Multigrate_s1.csv', index=False)