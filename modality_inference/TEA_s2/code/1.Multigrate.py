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

dir = "./imputation/input/TEA/"
n_batches = 4
rnas = []
atacs = []
adts = []

rna_input = []
atac_input = []
adt_input = []

rna_t = [1,2]
atac_t = [0,1]
adt_t =[2,3]

for batch in range(n_batches):
    
    if batch in rna_t:
        rds_data = pyreadr.read_r(os.path.join(dir, 'B' + str(batch + 2) + "/rna.rds"))
        counts_rna = np.array(rds_data[None]).T
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch + 1), counts_rna.shape[0])
        obs['modality'] = np.repeat('batch' + str(batch + 1), counts_rna.shape[0])
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
        rds_data = pyreadr.read_r(os.path.join(dir, 'B' + str(batch + 2) + "/atac.rds"))
        counts_atac = np.array(rds_data[None]).T
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch + 1), counts_atac.shape[0])
        obs['modality'] = np.repeat('batch' + str(batch + 1), counts_atac.shape[0])
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
        rds_data = pyreadr.read_r(os.path.join(dir, 'B' + str(batch + 2) + "/adt.rds"))
        counts_protein = np.array(rds_data[None]).T
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch + 1), counts_protein.shape[0])
        obs['modality'] = np.repeat('batch' + str(batch + 1), counts_protein.shape[0])
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

adata = mtg.data.organize_multiome_anndatas(
    adatas = [[rnas[0],rnas[1],rnas[2],rnas[3]], 
              [atacs[0],atacs[1],atacs[2],atacs[3]], 
              [adts[0],adts[1],adts[2],adts[3],]],
    layers = [[rna_input[0],rna_input[1],rna_input[2],rna_input[3]], 
              [atac_input[0],atac_input[1],atac_input[2],atac_input[3]], 
              [adt_input[0],adt_input[1],adt_input[2],adt_input[3]]],
)

del rnas
del atacs
del adts
del rds_data
del counts_protein
del counts_atac
del counts_rna
import gc
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
# atac
infer_res2 = model.impute(target_modality=1)
# adt
infer_res3 = model.impute(target_modality=2)


# B2 rna
B2_rna = infer_res1[0:6194,:]
if isinstance(B2_rna, np.ndarray):
    latent_df = pd.DataFrame(B2_rna)
else:
    
    latent_df = B2_rna

latent_df.to_csv('./imputation/TEA_s2/Multigrate/B2_rna.csv', index=False)

# B2 adt
B2_adt = infer_res3[0:6194,:]
if isinstance(B2_adt, np.ndarray):
    latent_df = pd.DataFrame(B2_adt)
else:
    
    latent_df = B2_adt

latent_df.to_csv('./imputation/TEA_s2/Multigrate/B2_adt.csv', index=False)

# B3 adt
B3_adt = infer_res3[6194:12575,:]
if isinstance(B3_adt, np.ndarray):
    latent_df = pd.DataFrame(B3_adt)
else:
    
    latent_df = B3_adt

latent_df.to_csv('./imputation/TEA_s2/Multigrate/B3_adt.csv', index=False)


# B4 atac
from scipy import sparse
from scipy import io
B4_atac = infer_res2[12575:18987,:]
csr_matrix = sparse.csr_matrix(B4_atac)
io.mmwrite("./imputation/TEA_s2/Multigrate/B4_atac.mtx", csr_matrix)

# B5 atac
B5_atac = infer_res2[18987:25517,:]
csr_matrix = sparse.csr_matrix(B5_atac)
io.mmwrite("./imputation/TEA_s2/Multigrate/B5_atac.mtx", csr_matrix)

# B5 rna
B5_rna = infer_res1[18987:25517,:]
if isinstance(B5_rna, np.ndarray):
    latent_df = pd.DataFrame(B5_rna)
else:
    
    latent_df = B5_rna

latent_df.to_csv('./imputation/TEA_s2/Multigrate/B5_rna.csv', index=False)
