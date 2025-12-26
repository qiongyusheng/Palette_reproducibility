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
mouses = []
humans = []

rna_input = []
mouse_input = []
human_input = []

rna_t = range(n_batches)
human_t = range(22)
mouse_t = range(22, 22 + 14)

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
        if batch in human_t:
            obs['modality'] = np.repeat('human', counts_rna.shape[0])
        else:
            obs['modality'] = np.repeat('mouse', counts_rna.shape[0])
        # np.repeat('rna', counts_rna.shape[0])
        rna_ad = ad.AnnData(counts_rna ,obs = obs)
        rna_ad.obs.index = rds_data[None].keys()
        rna_ad.layers["counts"] = rna_ad.X.copy()
        rna_in = "counts"
        
        # counts_rna = pyreadr.read_r(os.path.join(dir, 'rna' + str(batch + 1) + ".rds")).toarray().T
    else:
        rna_ad = None
        rna_in = None
    
    
    if batch in human_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/non.rds"))
        counts_rna = np.array(rds_data[None]).T
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat(path[batch], counts_rna.shape[0])
        obs['modality'] = np.repeat('human', counts_rna.shape[0])

        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_protein.shape[0])
        # np.repeat('adt', counts_protein.shape[0])
        human_ad = ad.AnnData(counts_rna ,obs = obs)
        human_ad.obs.index = rds_data[None].keys()
        sc.pp.normalize_total(human_ad)
        sc.pp.log1p(human_ad)
        human_ad.layers["counts"] = human_ad.X.copy()
        human_in = "counts"
        
    else:
        human_ad = None
        human_in = None
        
    if batch in mouse_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/non.rds"))
        counts_rna = np.array(rds_data[None]).T
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat(path[batch], counts_rna.shape[0])
        obs['modality'] = np.repeat('mouse', counts_rna.shape[0])

        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_protein.shape[0])
        # np.repeat('adt', counts_protein.shape[0])
        mouse_ad = ad.AnnData(counts_rna ,obs = obs)
        mouse_ad.obs.index = rds_data[None].keys()
        sc.pp.normalize_total(mouse_ad)
        sc.pp.log1p(mouse_ad)
        mouse_ad.layers["counts"] = mouse_ad.X.copy()
        mouse_in = "counts"
        
    else:
        mouse_ad = None
        mouse_in = None     
        
        
    rnas.append(rna_ad)
    rna_input.append(rna_in)
    
    humans.append(human_ad)
    human_input.append(human_in)
    
    mouses.append(mouse_ad)
    mouse_input.append(mouse_in)
    
adata = sca.models.organize_multiome_anndatas(
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
              rnas[33],rnas[34],rnas[35]], 
              [humans[0],humans[1],humans[2],
              humans[3],humans[4],humans[5],
              humans[6],humans[7],humans[8],
              humans[9],humans[10],humans[11],
              humans[12],humans[13],humans[14],
              humans[15],humans[16],humans[17],
              humans[18],humans[19],humans[20],
              humans[21],humans[22],humans[23],
              humans[24],humans[25],humans[26],
              humans[27],humans[28],humans[29],
              humans[30],humans[31],humans[32],
              humans[33],humans[34],humans[35]], 
              [mouses[0],mouses[1],mouses[2],
              mouses[3],mouses[4],mouses[5],
              mouses[6],mouses[7],mouses[8],
              mouses[9],mouses[10],mouses[11],
              mouses[12],mouses[13],mouses[14],
              mouses[15],mouses[16],mouses[17],
              mouses[18],mouses[19],mouses[20],
              mouses[21],mouses[22],mouses[23],
              mouses[24],mouses[25],mouses[26],
              mouses[27],mouses[28],mouses[29],
              mouses[30],mouses[31],mouses[32],
              mouses[33],mouses[34],mouses[35]]],
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
              rna_input[33],rna_input[34],rna_input[35]], 
              
              [human_input[0],human_input[1],human_input[2],
              human_input[3],human_input[4],human_input[5],
              human_input[6],human_input[7],human_input[8],
              human_input[9],human_input[10],human_input[11],
              human_input[12],human_input[13],human_input[14],
              human_input[15],human_input[16],human_input[17],
              human_input[18],human_input[19],human_input[20],
              human_input[21],human_input[22],human_input[23],
              human_input[24],human_input[25],human_input[26],
              human_input[27],human_input[28],human_input[29],
              human_input[30],human_input[31],human_input[32],
              human_input[33],human_input[34],human_input[35]],
              
              [mouse_input[0],mouse_input[1],mouse_input[2],
              mouse_input[3],mouse_input[4],mouse_input[5],
              mouse_input[6],mouse_input[7],mouse_input[8],
              mouse_input[9],mouse_input[10],mouse_input[11],
              mouse_input[12],mouse_input[13],mouse_input[14],
              mouse_input[15],mouse_input[16],mouse_input[17],
              mouse_input[18],mouse_input[19],mouse_input[20],
              mouse_input[21],mouse_input[22],mouse_input[23],
              mouse_input[24],mouse_input[25],mouse_input[26],
              mouse_input[27],mouse_input[28],mouse_input[29],
              mouse_input[30],mouse_input[31],mouse_input[32],
              mouse_input[33],mouse_input[34],mouse_input[35]]],
)

del rnas
del humans
del mouses
del rds_data
del counts_rna
import gc
gc.collect()

sca.models.MultiVAE.setup_anndata(
    adata,
    categorical_covariate_keys=['modality', 'samplename'],
    rna_indices_end=4000,
)

model = sca.models.MultiVAE(
    adata,
    losses=['nb', 'mse', 'mse'],
    loss_coefs={'kl': 1e-1,
               'integ': 3000,
               },
    integrate_on='modality',
    mmd='marginal',z_dim = 20,
)
model.train()
model.get_latent_representation()
latent_data = adata.obsm['latent']
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    latent_df = latent_data

latent_df.to_csv('./Multigrate_non.csv', 
                 index=False)

