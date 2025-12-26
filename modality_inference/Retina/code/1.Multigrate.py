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
import scipy.io as sio
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


dir = "./imputation/input/Retina/"
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
    
        
    rnas.append(rna_ad)
    rna_input.append(rna_in)
    
    atacs.append(atac_ad)
    atac_input.append(atac_in)
    
import multigrate as mtg

del rds_data
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

# rna
infer_res1 = model.impute(target_modality=0)
# atac
# infer_res2 = model.impute(target_modality=1)

# LGS3OS rna
LGS3OS_rna = infer_res1[35533:40546,:]

if isinstance(LGS3OS_rna, np.ndarray):
    latent_df = pd.DataFrame(LGS3OS_rna)
else:
    
    latent_df = LGS3OS_rna


latent_df.to_csv('./imputation/Retina/Multigrate/LGS3OS_rna.csv', 
                 index=False)
del LGS3OS_rna
del latent_df

# LVG1OD rna
LVG1OD_rna = infer_res1[40546:45583,:]

if isinstance(LVG1OD_rna, np.ndarray):
    latent_df = pd.DataFrame(LVG1OD_rna)
else:
    
    latent_df = LVG1OD_rna


latent_df.to_csv('./imputation/Retina/Multigrate/LVG1OD_rna.csv', 
                 index=False)
del LVG1OD_rna
del latent_df

# LVG1OS rna
LVG1OS_rna = infer_res1[45583:50312,:]

if isinstance(LVG1OS_rna, np.ndarray):
    latent_df = pd.DataFrame(LVG1OS_rna)
else:
    
    latent_df = LVG1OS_rna


latent_df.to_csv('./imputation/Retina/Multigrate/LVG1OS_rna.csv', 
                 index=False)

del LVG1OS_rna
del latent_df
del infer_res1

# atac
infer_res2 = model.impute(target_modality=1)
from scipy import sparse
from scipy import io

# LGS2OD atac
LGS2OD_atac = infer_res2[17637:23871,:]
csr_matrix = sparse.csr_matrix(LGS2OD_atac)
io.mmwrite("./imputation/Retina/Multigrate/LGS2OD_atac.mtx", csr_matrix)
del LGS2OD_atac
del csr_matrix

# LGS2OS atac
LGS2OS_atac = infer_res2[23871:29826,:]
csr_matrix = sparse.csr_matrix(LGS2OS_atac)
io.mmwrite("./imputation/Retina/Multigrate/LGS2OS_atac.mtx", csr_matrix)

del LGS2OS_atac
del csr_matrix

# LGS3OD atac
LGS3OD_atac = infer_res2[29826:35533,:]
csr_matrix = sparse.csr_matrix(LGS3OD_atac)
io.mmwrite("./imputation/Retina/Multigrate/LGS3OD_atac.mtx", csr_matrix)

del LGS3OD_atac 
del csr_matrix