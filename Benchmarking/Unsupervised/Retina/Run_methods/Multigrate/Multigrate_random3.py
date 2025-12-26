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

rna_t = [0,1,2,3,4]
atac_t = [0,1,5,6,7]

path = ['LGS1OD',
        'LGS1OS',
        'LGS2OD',
        'LGS2OS',
        'LGS3OD',
        'LGS3OS',
        'LVG1OD',
        'LVG1OS']

atac_cell_name = []

for batch in range(n_batches):
    
    if batch in rna_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/rna.rds"))
        counts_rna = np.array(rds_data[None]).T
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch), counts_rna.shape[0])
        if batch in [0,1]:
            obs['modality'] = np.repeat('Multiome', counts_rna.shape[0])
        else:
            obs['modality'] = np.repeat('RNA', counts_rna.shape[0])

        # np.repeat('rna', counts_rna.shape[0])
        rna_ad = ad.AnnData(counts_rna ,obs = obs)
        rna_ad.obs.index = rds_data[None].keys()
        rna_ad.layers["counts"] = rna_ad.X.copy()
        rna_in = "counts"
        

    else:
        rna_ad = None
        rna_in = None
    
    if batch in atac_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/atac.rds"))
        atac_ad = np.array(rds_data[None]).T
        atac_cell_name.append(rds_data[None].keys())
        if batch in [0]:
            peak_name = rds_data[None].index
        
        atac_in = "log-norm"
        
    else:
        atac_ad = None
        atac_in = None
    
        
    rnas.append(rna_ad)
    rna_input.append(rna_in)
    
    atacs.append(atac_ad)
    atac_input.append(atac_in)
    
tmp = np.concatenate((atac_ad[0], atac_ad[1], atac_ad[5],atac_ad[6],atac_ad[7]), axis=0)
adata = ad.AnnData(tmp,dtype=np.float32)
adata.obs_name = np.concatenate(atac_cell_name,axis = 0)
adata.var_name = peak_name
adata.obs['samplename'] = ["batch0"] * atac_ad[0].shape[0] + ["batch1"] * atac_ad[1].shape[0] + ["batch5"] * atac_ad[5].shape[0] + ["batch6"] * atac_ad[6].shape[0] + ["batch7"] * atac_ad[7].shape[0]
adata.obs['modality'] = ["Multiome"] * atac_ad[0].shape[0] + ["Multiome"] * atac_ad[1].shape[0] + ["ATAC"] * atac_ad[5].shape[0] + ["ATAC"] * atac_ad[6].shape[0] + ["ATAC"] * atac_ad[7].shape[0]
sc.pp.highly_variable_genes(adata, n_top_genes=50000, batch_key='samplename')
adata_hvg = adata[:, adata.var.highly_variable].copy()
adata_hvg.layers['log-norm'] = adata_hvg.X.copy()
adata_hvg0 = adata_hvg[adata.obs['samplename'] == 'batch0',:]
adata_hvg1 = adata_hvg[adata.obs['samplename'] == 'batch1',:]
adata_hvg5 = adata_hvg[adata.obs['samplename'] == 'batch5',:]
adata_hvg6 = adata_hvg[adata.obs['samplename'] == 'batch6',:]
adata_hvg7 = adata_hvg[adata.obs['samplename'] == 'batch7',:]
atacs[0] = adata_hvg0
atacs[1] = adata_hvg1
atacs[5] = adata_hvg5
atacs[6] = adata_hvg6
atacs[7] = adata_hvg7

import multigrate as mtg

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
latent_df.to_csv('./Benchmarking/Unsupervised/Retina/output/Multigrate.csv', index=False)