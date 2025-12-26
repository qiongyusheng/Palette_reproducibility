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

dir = "./imputation/input/BMMC_s/"
n_batches = 10
rnas = []
atacs = []
adts = []

rna_input = []
atac_input = []
adt_input = []

rna_t = [0,2,4,5,6,8]
atac_t = [5,7,9]
adt_t =[1,3,4]

path = ['CITE1rna','CITE1adt','CITE2rna','CITE2adt','CITE3',
        'Multiome1','Multiome2rna','Multiome2atac','Multiome3rna','Multiome3atac']

for batch in range(n_batches):
    
    if batch in rna_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/rna.rds"))
        counts_rna = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch), counts_rna.shape[0])
        if batch in [4]:
            obs['modality'] = np.repeat('CITE', counts_rna.shape[0])
        else:
            if batch in [5]:
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
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/atac.rds"))
        peak_name = np.array(rds_data[None].index, dtype=str)
        counts_atac = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch), counts_atac.shape[0])
        # obs['modality'] = np.repeat('Multiome', counts_atac.shape[0])
        if batch in [5]:
            obs['modality'] = np.repeat('Multiome', counts_atac.shape[0])
        else:
            obs['modality'] = np.repeat('ATAC', counts_atac.shape[0])
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
        counts_protein = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch), counts_protein.shape[0])
        if batch in [4]:
            obs['modality'] = np.repeat('CITE', counts_protein.shape[0])
        else:
            obs['modality'] = np.repeat('ADT', counts_protein.shape[0])
        # obs['modality'] = np.repeat('CITE', counts_protein.shape[0])
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

# rna
infer_res1 = model.impute(target_modality=0)
rna_adata = ad.AnnData(infer_res1)
rna_adata.obs['samplename'] = np.array(adata.obs['samplename'], dtype=str)
del infer_res1

# cite1 rna
adata_cite1 = rna_adata[rna_adata.obs['samplename'] == 'batch1',:]
cite1_rna = adata_cite1.X

if isinstance(cite1_rna, np.ndarray):
    latent_df = pd.DataFrame(cite1_rna)
else:
    
    latent_df = cite1_rna
    
latent_df.to_csv('./imputation/BMMC_s/Multigrate/cite1rna.csv',
                 index=False)  
del adata_cite1
del cite1_rna
del latent_df


# cite2 rna
adata_cite2 = rna_adata[rna_adata.obs['samplename'] == 'batch3',:]
cite2_rna = adata_cite2.X

if isinstance(cite2_rna, np.ndarray):
    latent_df = pd.DataFrame(cite2_rna)
else:
    
    latent_df = cite2_rna
    
latent_df.to_csv('./imputation/BMMC_s/Multigrate/cite2rna.csv',
                 index=False)  
del adata_cite2
del cite2_rna
del latent_df

# mul2 rna
adata_mul2 = rna_adata[rna_adata.obs['samplename'] == 'batch7',:]
mul2_rna = adata_mul2.X

if isinstance(mul2_rna, np.ndarray):
    latent_df = pd.DataFrame(mul2_rna)
else:
    
    latent_df = mul2_rna
    
latent_df.to_csv('./imputation/BMMC_s/Multigrate/mul2rna.csv',
                 index=False)  
del adata_mul2
del mul2_rna
del latent_df

# mul3 rna
adata_mul3 = rna_adata[rna_adata.obs['samplename'] == 'batch9',:]
mul3_rna = adata_mul3.X

if isinstance(mul3_rna, np.ndarray):
    latent_df = pd.DataFrame(mul3_rna)
else:
    
    latent_df = mul3_rna
    
latent_df.to_csv('./imputation/BMMC_s/Multigrate/mul3rna.csv',
                 index=False)  
del adata_mul3
del mul3_rna
del latent_df

del rna_adata

# adt
infer_res3 = model.impute(target_modality=2)
adt_adata = ad.AnnData(infer_res3)
adt_adata.obs['samplename'] = np.array(adata.obs['samplename'], dtype=str)
del infer_res3

# cite1 adt
adata_cite1 = adt_adata[adt_adata.obs['samplename'] == 'batch0',:]
cite1_adt = adata_cite1.X

if isinstance(cite1_adt, np.ndarray):
    latent_df = pd.DataFrame(cite1_adt)
else:
    
    latent_df = cite1_adt
    
latent_df.to_csv('./imputation/BMMC_s/Multigrate/cite1adt.csv',
                 index=False)  
del adata_cite1
del cite1_adt
del latent_df

# cite2 adt
adata_cite2 = adt_adata[adt_adata.obs['samplename'] == 'batch2',:]
cite2_adt = adata_cite2.X

if isinstance(cite2_adt, np.ndarray):
    latent_df = pd.DataFrame(cite2_adt)
else:
    
    latent_df = cite2_adt
    
latent_df.to_csv('./imputation/BMMC_s/Multigrate/cite2adt.csv',
                 index=False)  
del adata_cite2
del cite2_adt
del latent_df

del adt_adata

# atac
infer_res2 = model.impute(target_modality=1)
infer_res2_csr = csr_matrix(infer_res2)
del infer_res2
atac_adata = ad.AnnData(infer_res2_csr)
atac_adata.obs['samplename'] = np.array(adata.obs['samplename'], dtype=str)
del infer_res2_csr

# mul2 atac
adata_mul2 = atac_adata[atac_adata.obs['samplename'] == 'batch6',:]
mul2_atac = adata_mul2.X

io.mmwrite("./imputation/BMMC_s/Multigrate/mul2atac.mtx", 
           mul2_atac)
del adata_mul2
del mul2_atac

# mul3 atac
adata_mul3 = atac_adata[atac_adata.obs['samplename'] == 'batch8',:]
mul3_atac = adata_mul3.X

io.mmwrite("./imputation/BMMC_s/Multigrate/mul3atac.mtx", 
           mul3_atac)
del adata_mul3
del mul3_atac