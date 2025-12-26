# # 1.Reference model train
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

import warnings
warnings.filterwarnings("ignore")

dir = "./BMMC_CITE_Multiome/"
n_batches = 6
rnas = []
atacs = []
adts = []

rna_input = []
atac_input = []
adt_input = []

rna_t = [0,1,2,3,4,5]
atac_t = [3,4,5]
adt_t =[0,1,2]

path = ['CITE/s1d1_cite','CITE/s1d2_cite','CITE/s1d3_cite',
        'Multiome/s1d1_multi','Multiome/s1d2_multi','Multiome/s1d3_multi']

for batch in range(n_batches):
    
    if batch in rna_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/rna.rds"))
        counts_rna = np.array(rds_data[None]).T
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
        rna_in = None
        
        # counts_rna = pyreadr.read_r(os.path.join(dir, 'rna' + str(batch + 1) + ".rds")).toarray().T
    else:
        rna_ad = None
        rna_in = None
    
    if batch in atac_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/atac.rds"))
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
        atac_in = None
        
    else:
        atac_ad = None
        atac_in = None
    
    if batch in adt_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/adt.rds"))
        counts_protein = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch), counts_protein.shape[0])
        obs['modality'] = np.repeat('CITE', counts_protein.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_protein.shape[0])
        # np.repeat('adt', counts_protein.shape[0])
        adt_ad = ad.AnnData(counts_protein ,obs = obs)
        adt_ad.obs.index = rds_data[None].keys()
        muon.prot.pp.clr(adt_ad)
        adt_in = None
        
    else:
        adt_ad = None
        adt_in = None
        
    rnas.append(rna_ad)
    rna_input.append(rna_in)
    
    atacs.append(atac_ad)
    atac_input.append(atac_in)
    
    adts.append(adt_ad)
    adt_input.append(adt_in)

del rds_data
del counts_protein
del counts_atac
del counts_rna
import gc
gc.collect()

import multigrate as mtg

adata = mtg.data.organize_multiome_anndatas(
    adatas = [[rnas[0],rnas[1],rnas[2],rnas[3],rnas[4],rnas[5]], 
              [atacs[0],atacs[1],atacs[2],atacs[3],atacs[4],atacs[5]], 
              [adts[0],adts[1],adts[2],adts[3],adts[4],adts[5]]],
    layers = [[rna_input[0],rna_input[1],rna_input[2],rna_input[3],rna_input[4],rna_input[5]], 
              [atac_input[0],atac_input[1],atac_input[2],atac_input[3],atac_input[4],atac_input[5]], 
              [adt_input[0],adt_input[1],adt_input[2],adt_input[3],adt_input[4],adt_input[5]]],
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

models = mtg.model.MultiVAE(
    adata,
    losses=['nb', 'mse', 'mse'],
    loss_coefs={'kl': 1e-3,
               'integ': 3000,
               },
    integrate_on='modality',
    mmd='marginal',
    z_dim = 20,
)

models.train()
models.get_latent_representation()


path = ['CITE/s2d1_cite','CITE/s2d4_cite','CITE/s2d5_cite',
       'Multiome/s2d1_multi','Multiome/s2d4_multi','Multiome/s2d5_multi']

rds_rna = pyreadr.read_r(os.path.join(dir, path[0] + "/rna.rds"))

rds_adt = pyreadr.read_r(os.path.join(dir, path[0] + "/adt.rds"))

rds_atac = pyreadr.read_r(os.path.join(dir, path[3] + "/atac.rds"))


# # ADT transfer

n_batches = 3
adts = []
rnas = []
atacs = []

adt_input = []
rna_input = []
atac_input = []

adt_t =[0,1,2]
rna_t =[0,1,2]
atac_t =[0,1,2]

path = ['CITE/s2d1_cite','CITE/s2d4_cite','CITE/s2d5_cite']

for batch in range(n_batches):
    
    if batch in adt_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/adt.rds"))
        counts_protein = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_protein.shape[0])
        obs['modality'] = np.repeat('ADT', counts_protein.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_protein.shape[0])
        # np.repeat('adt', counts_protein.shape[0])
        adt_ad = ad.AnnData(counts_protein ,obs = obs)
        adt_ad.obs.index = rds_data[None].keys()
        muon.prot.pp.clr(adt_ad)
        # adt_ad.layers['clr'] = adt_ad.X.copy()
        # adt_in = "clr"
        adt_in = None
        
    else:
        adt_ad = None
        adt_in = None
    
    if batch in rna_t:
        # rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/rna.rds"))
        counts_rna = csr_matrix((counts_protein.shape[0], rds_rna[None].shape[0]))
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_rna.shape[0])
        obs['modality'] = np.repeat('ADT', counts_rna.shape[0])
        # np.repeat('rna', counts_rna.shape[0])
        rna_ad = ad.AnnData(counts_rna ,obs = obs)
        rna_ad.obs.index = adt_ad.obs.index
        # rna_ad.layers["counts"] = rna_ad.X.copy()
        # rna_in = "counts"
        rna_in = None
        
        # counts_rna = pyreadr.read_r(os.path.join(dir, 'rna' + str(batch + 1) + ".rds")).toarray().T
    else:
        rna_ad = None
        rna_in = None
    
    if batch in atac_t:
        # rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/atac.rds"))
        counts_atac = csr_matrix((counts_protein.shape[0], rds_atac[None].shape[0]))

        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_atac.shape[0])
        obs['modality'] = np.repeat('ADT', counts_atac.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_atac.shape[0])
        # np.repeat('atac', counts_atac.shape[0])
        atac_ad = ad.AnnData(counts_atac ,obs = obs)
        atac_ad.obs.index = adt_ad.obs.index
        atac_in = None
        
    else:
        atac_ad = None
        atac_in = None
    
        
    rnas.append(rna_ad)
    rna_input.append(rna_in)
    
    atacs.append(atac_ad)
    atac_input.append(atac_in)
    
    adts.append(adt_ad)
    adt_input.append(adt_in)

del rds_data
del counts_protein
del counts_atac
del counts_rna
import gc
gc.collect()

query = mtg.data.organize_multiome_anndatas(
    adatas = [[rnas[0],rnas[1],rnas[2]], 
              [atacs[0],atacs[1],atacs[2]], 
              [adts[0],adts[1],adts[2]]],
    layers = [[rna_input[0],rna_input[1],rna_input[2]], 
              [atac_input[0],atac_input[1],atac_input[2]], 
              [adt_input[0],adt_input[1],adt_input[2]]],
)
del rnas
del atacs
del adts
gc.collect()

new_model = mtg.model.MultiVAE.load_query_data(query,models)

new_model.train(weight_decay=0)
new_model.get_latent_representation()

adata_both = ad.concat([adata, query])
latent_data = adata_both.obsm['latent']

if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:

    latent_df = latent_data


latent_df.to_csv('/home/server/sqy/data/mosaic/result/mapping/adt/Multigrate.csv', index=False)


del latent_df
del latent_data
del adata_both
del query
del new_model


# # RNA transfer

n_batches = 6
rnas = []
adts = []
atacs = []

rna_input = []
adt_input = []
atac_input = []

rna_t = [0,1,2,3,4,5]
adt_t = [0,1,2,3,4,5]
atac_t = [0,1,2,3,4,5]

path = ['CITE/s2d1_cite','CITE/s2d4_cite','CITE/s2d5_cite',
       'Multiome/s2d1_multi','Multiome/s2d4_multi','Multiome/s2d5_multi']


for batch in range(n_batches):
    
    if batch in rna_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/rna.rds"))
        counts_rna = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_rna.shape[0])
        obs['modality'] = np.repeat('RNA', counts_rna.shape[0])
        rna_ad = ad.AnnData(counts_rna ,obs = obs)
        rna_ad.obs.index = rds_data[None].keys()
        rna_in = None
        
        # counts_rna = pyreadr.read_r(os.path.join(dir, 'rna' + str(batch + 1) + ".rds")).toarray().T
    else:
        rna_ad = None
        rna_in = None
    
    if batch in adt_t:
        counts_protein = csr_matrix((counts_rna.shape[0], rds_adt[None].shape[0]))
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_protein.shape[0])
        obs['modality'] = np.repeat('RNA', counts_protein.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_protein.shape[0])
        # np.repeat('adt', counts_protein.shape[0])
        adt_ad = ad.AnnData(counts_protein ,obs = obs)
        adt_ad.obs.index = rna_ad.obs.index
        muon.prot.pp.clr(adt_ad)
        # adt_ad.layers['clr'] = adt_ad.X.copy()
        # adt_in = "clr"
        adt_in = None
        
    else:
        adt_ad = None
        adt_in = None
        
    if batch in atac_t:
        # rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/atac.rds"))
        counts_atac = csr_matrix((counts_rna.shape[0], rds_atac[None].shape[0]))

        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_atac.shape[0])
        obs['modality'] = np.repeat('RNA', counts_atac.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_atac.shape[0])
        # np.repeat('atac', counts_atac.shape[0])
        atac_ad = ad.AnnData(counts_atac ,obs = obs)
        atac_ad.obs.index = rna_ad.obs.index
        atac_in = None
        
    else:
        atac_ad = None
        atac_in = None  
        
        
    rnas.append(rna_ad)
    rna_input.append(rna_in)
    
    atacs.append(atac_ad)
    atac_input.append(atac_in)
    
    adts.append(adt_ad)
    adt_input.append(adt_in)

del rds_data
del counts_protein
del counts_atac
del counts_rna
import gc
gc.collect()

query = mtg.data.organize_multiome_anndatas(
    adatas = [[rnas[0],rnas[1],rnas[2],rnas[3],rnas[4],rnas[5]], 
              [atacs[0],atacs[1],atacs[2],atacs[3],atacs[4],atacs[5]], 
              [adts[0],adts[1],adts[2],adts[3],adts[4],adts[5]]],
    layers = [[rna_input[0],rna_input[1],rna_input[2],rna_input[3],rna_input[4],rna_input[5]], 
              [atac_input[0],atac_input[1],atac_input[2],atac_input[3],atac_input[4],atac_input[5]], 
              [adt_input[0],adt_input[1],adt_input[2],adt_input[3],adt_input[4],adt_input[5]]],
)
del rnas
del atacs
del adts
gc.collect()


new_model = mtg.model.MultiVAE.load_query_data(query,models)
new_model.train(weight_decay=0)
new_model.get_latent_representation()

adata_both = ad.concat([adata, query])
latent_data = adata_both.obsm['latent']

if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:

    latent_df = latent_data

latent_df.to_csv('/home/server/sqy/data/mosaic/result/mapping/rna/Multigrate.csv', index=False)


del latent_df
del latent_data
del adata_both
del query
del new_model


# # RNA + ATAC transfer

n_batches = 6
rnas = []
atacs = []
adts = []

rna_input = []
adt_input = []
atac_input = []

rna_t = [0,1,2]
atac_t = [3,4,5]
adt_t = [0,1,2]

path = ['CITE/s2d1_cite','CITE/s2d4_cite','CITE/s2d5_cite',
       'Multiome/s2d1_multi','Multiome/s2d4_multi','Multiome/s2d5_multi']


for batch in range(n_batches):
    
    if batch in rna_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/rna.rds"))
        counts_rna = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_rna.shape[0])
        obs['modality'] = np.repeat('RNA', counts_rna.shape[0])
        # np.repeat('rna', counts_rna.shape[0])
        rna_ad = ad.AnnData(counts_rna ,obs = obs)
        rna_ad.obs.index = rds_data[None].keys()
        rna_in = None
        
        # counts_rna = pyreadr.read_r(os.path.join(dir, 'rna' + str(batch + 1) + ".rds")).toarray().T
    else:
        rna_ad = None
        rna_in = None
    
    if batch in adt_t:
        counts_protein = csr_matrix((counts_rna.shape[0],rds_adt[None].shape[0]))
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_protein.shape[0])
        obs['modality'] = np.repeat('RNA', counts_protein.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_protein.shape[0])
        # np.repeat('adt', counts_protein.shape[0])
        adt_ad = ad.AnnData(counts_protein ,obs = obs)
        adt_ad.obs.index = rna_ad.obs.index
        adt_in = None
        
    else:
        adt_ad = None
        adt_in = None
    
    if batch in atac_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/atac.rds"))
        counts_atac = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_atac.shape[0])
        obs['modality'] = np.repeat('ATAC', counts_atac.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_atac.shape[0])
        # np.repeat('atac', counts_atac.shape[0])
        atac_ad = ad.AnnData(counts_atac ,obs = obs)
        atac_ad.obs.index = rds_data[None].keys()
        sc.pp.normalize_total(atac_ad, target_sum=1e4)
        sc.pp.log1p(atac_ad)
        atac_in = None
        
    else:
        atac_ad = None
        atac_in = None
    
        
    rnas.append(rna_ad)
    rna_input.append(rna_in)
    
    atacs.append(atac_ad)
    atac_input.append(atac_in)
    
    adts.append(adt_ad)
    adt_input.append(adt_in)


del rds_data
del counts_protein
del counts_atac
del counts_rna
import gc
gc.collect()

query = mtg.data.organize_multiome_anndatas(
    adatas = [[rnas[0],rnas[1],rnas[2],rnas[3],rnas[4],rnas[5]], 
              [atacs[0],atacs[1],atacs[2],atacs[3],atacs[4],atacs[5]], 
              [adts[0],adts[1],adts[2],adts[3],adts[4],adts[5]]],
    layers = [[rna_input[0],rna_input[1],rna_input[2],rna_input[3],rna_input[4],rna_input[5]], 
              [atac_input[0],atac_input[1],atac_input[2],atac_input[3],atac_input[4],atac_input[5]], 
              [adt_input[0],adt_input[1],adt_input[2],adt_input[3],adt_input[4],adt_input[5]]],
)


del rnas
del atacs
del adts
gc.collect()


new_model = mtg.model.MultiVAE.load_query_data(query,models)
new_model.train(weight_decay=0)
new_model.get_latent_representation()


adata_both = ad.concat([adata, query])
latent_data = adata_both.obsm['latent']

if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    latent_df = latent_data


latent_df.to_csv('/home/server/sqy/data/mosaic/result/mapping/rna_atac/Multigrate.csv', index=False)


# # RNA + ADT Transfer

n_batches = 6
rnas = []
adts = []
atacs = []

rna_input = []
atac_input = []
adt_input = []

rna_t = [3,4,5]
atac_t = [3,4,5]
adt_t =[0,1,2]

path = ['CITE/s2d1_cite','CITE/s2d4_cite','CITE/s2d5_cite',
       'Multiome/s2d1_multi','Multiome/s2d4_multi','Multiome/s2d5_multi']


for batch in range(n_batches):
    
    if batch in rna_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/rna.rds"))
        counts_rna = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_rna.shape[0])
        obs['modality'] = np.repeat('RNA', counts_rna.shape[0])
        # np.repeat('rna', counts_rna.shape[0])
        rna_ad = ad.AnnData(counts_rna ,obs = obs)
        rna_ad.obs.index = rds_data[None].keys()
        rna_in = None
        
        # counts_rna = pyreadr.read_r(os.path.join(dir, 'rna' + str(batch + 1) + ".rds")).toarray().T
    else:
        rna_ad = None
        rna_in = None
    
    
    if batch in adt_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/adt.rds"))
        counts_protein = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_protein.shape[0])
        obs['modality'] = np.repeat('ADT', counts_protein.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_protein.shape[0])
        # np.repeat('adt', counts_protein.shape[0])
        adt_ad = ad.AnnData(counts_protein ,obs = obs)
        adt_ad.obs.index = rds_data[None].keys()
        muon.prot.pp.clr(adt_ad)
        rna_in = None
        
    else:
        adt_ad = None
        adt_in = None
        
    if batch in atac_t:
        counts_atac = csr_matrix((counts_rna.shape[0],rds_atac[None].shape[0]))
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_atac.shape[0])
        obs['modality'] = np.repeat('RNA', counts_atac.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_atac.shape[0])
        # np.repeat('atac', counts_atac.shape[0])
        atac_ad = ad.AnnData(counts_atac ,obs = obs)
        atac_ad.obs.index = rna_ad.obs.index
        atac_in = None
        
    else:
        atac_ad = None
        atac_in = None

        
    rnas.append(rna_ad)
    rna_input.append(rna_in)
    
    atacs.append(atac_ad)
    atac_input.append(atac_in)
    
    adts.append(adt_ad)
    adt_input.append(adt_in)


del rds_data
del counts_protein
del counts_atac
del counts_rna
import gc
gc.collect()


query = mtg.data.organize_multiome_anndatas(
    adatas = [[rnas[0],rnas[1],rnas[2],rnas[3],rnas[4],rnas[5]], 
              [atacs[0],atacs[1],atacs[2],atacs[3],atacs[4],atacs[5]], 
              [adts[0],adts[1],adts[2],adts[3],adts[4],adts[5]]],
    layers = [[rna_input[0],rna_input[1],rna_input[2],rna_input[3],rna_input[4],rna_input[5]], 
              [atac_input[0],atac_input[1],atac_input[2],atac_input[3],atac_input[4],atac_input[5]], 
              [adt_input[0],adt_input[1],adt_input[2],adt_input[3],adt_input[4],adt_input[5]]],
)
del rnas
del atacs
del adts
gc.collect()


new_model = mtg.model.MultiVAE.load_query_data(query,models)
new_model.train(weight_decay=0)
new_model.get_latent_representation()


adata_both = ad.concat([adata, query])
latent_data = adata_both.obsm['latent']

if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:

    latent_df = latent_data


latent_df.to_csv('/home/server/sqy/data/mosaic/result/mapping/rna_adt/Multigrate.csv', index=False)


# # RNA + ADT + ATAC Transfer

n_batches = 6
rnas = []
atacs = []
adts = []

rna_input = []
atac_input = []
adt_input = []

rna_t = [2,3]
atac_t = [4,5]
adt_t =[0,1]

path = ['CITE/s2d1_cite','CITE/s2d4_cite','CITE/s2d5_cite',
       'Multiome/s2d1_multi','Multiome/s2d4_multi','Multiome/s2d5_multi']


for batch in range(n_batches):
    
    if batch in rna_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/rna.rds"))
        counts_rna = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_rna.shape[0])
        obs['modality'] = np.repeat('RNA', counts_rna.shape[0])
        # np.repeat('rna', counts_rna.shape[0])
        rna_ad = ad.AnnData(counts_rna ,obs = obs)
        rna_ad.obs.index = rds_data[None].keys()
        rna_in = None
        
        # counts_rna = pyreadr.read_r(os.path.join(dir, 'rna' + str(batch + 1) + ".rds")).toarray().T
    else:
        rna_ad = None
        rna_in = None
    
    if batch in atac_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/atac.rds"))
        counts_atac = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_atac.shape[0])
        obs['modality'] = np.repeat('ATAC', counts_atac.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_atac.shape[0])
        # np.repeat('atac', counts_atac.shape[0])
        atac_ad = ad.AnnData(counts_atac ,obs = obs)
        atac_ad.obs.index = rds_data[None].keys()
        sc.pp.normalize_total(atac_ad, target_sum=1e4)
        sc.pp.log1p(atac_ad)
        atac_in = None
        
    else:
        atac_ad = None
        atac_in = None
    
    if batch in adt_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/adt.rds"))
        counts_protein = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_protein.shape[0])
        obs['modality'] = np.repeat('ADT', counts_protein.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_protein.shape[0])
        # np.repeat('adt', counts_protein.shape[0])
        adt_ad = ad.AnnData(counts_protein ,obs = obs)
        adt_ad.obs.index = rds_data[None].keys()
        muon.prot.pp.clr(adt_ad)
        adt_in = None
        
    else:
        adt_ad = None
        adt_in = None
        
    rnas.append(rna_ad)
    rna_input.append(rna_in)
    
    atacs.append(atac_ad)
    atac_input.append(atac_in)
    
    adts.append(adt_ad)
    adt_input.append(adt_in)


del rds_data
del counts_protein
del counts_atac
del counts_rna
import gc
gc.collect()

query = mtg.data.organize_multiome_anndatas(
    adatas = [[rnas[0],rnas[1],rnas[2],rnas[3],rnas[4],rnas[5]], 
              [atacs[0],atacs[1],atacs[2],atacs[3],atacs[4],atacs[5]], 
              [adts[0],adts[1],adts[2],adts[3],adts[4],adts[5]]],
    layers = [[rna_input[0],rna_input[1],rna_input[2],rna_input[3],rna_input[4],rna_input[5]], 
              [atac_input[0],atac_input[1],atac_input[2],atac_input[3],atac_input[4],atac_input[5]], 
              [adt_input[0],adt_input[1],adt_input[2],adt_input[3],adt_input[4],adt_input[5]]],
)


del rnas
del atacs
del adts
gc.collect()


new_model = mtg.model.MultiVAE.load_query_data(query,models)
new_model.train(weight_decay=0)
new_model.get_latent_representation()


adata_both = ad.concat([adata, query])
latent_data = adata_both.obsm['latent']

if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:

    latent_df = latent_data


latent_df.to_csv('/home/server/sqy/data/mosaic/result/mapping/rna_adt_atac/Multigrate.csv', index=False)


# # Multiome Transfer

n_batches = 3
rnas = []
atacs = []
adts = []

rna_input = []
atac_input = []
adt_input = []

rna_t = [0,1,2]
atac_t = [0,1,2]
adt_t = [0,1,2]

path = [
       'Multiome/s2d1_multi','Multiome/s2d4_multi','Multiome/s2d5_multi']


for batch in range(n_batches):
    
    if batch in rna_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/rna.rds"))
        counts_rna = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_rna.shape[0])
        obs['modality'] = np.repeat('Multiome', counts_rna.shape[0])
        # np.repeat('rna', counts_rna.shape[0])
        rna_ad = ad.AnnData(counts_rna ,obs = obs)
        rna_ad.obs.index = rds_data[None].keys()
        rna_in = None
        
        # counts_rna = pyreadr.read_r(os.path.join(dir, 'rna' + str(batch + 1) + ".rds")).toarray().T
    else:
        rna_ad = None
        rna_in = None
    
    if batch in atac_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/atac.rds"))
        counts_atac = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_atac.shape[0])
        obs['modality'] = np.repeat('Multiome', counts_atac.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_atac.shape[0])
        # np.repeat('atac', counts_atac.shape[0])
        atac_ad = ad.AnnData(counts_atac ,obs = obs)
        atac_ad.obs.index = rds_data[None].keys()
        sc.pp.normalize_total(atac_ad, target_sum=1e4)
        sc.pp.log1p(atac_ad)
        atac_in = None
        
    else:
        atac_ad = None
        atac_in = None

    if batch in adt_t:
        counts_protein = csr_matrix((counts_atac.shape[0],rds_adt[None].shape[0]))
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_protein.shape[0])
        obs['modality'] = np.repeat('Multiome', counts_protein.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_protein.shape[0])
        # np.repeat('adt', counts_protein.shape[0])
        adt_ad = ad.AnnData(counts_protein ,obs = obs)
        adt_ad.obs.index = atac_ad.obs.index
        adt_in = None
        
    else:
        adt_ad = None
        adt_in = None
    
        
    rnas.append(rna_ad)
    rna_input.append(rna_in)
    
    atacs.append(atac_ad)
    atac_input.append(atac_in)
    
    adts.append(adt_ad)
    adt_input.append(adt_in)


del rds_data
del counts_protein
del counts_atac
del counts_rna
import gc
gc.collect()


query = mtg.data.organize_multiome_anndatas(
    adatas = [[rnas[0],rnas[1],rnas[2]], 
              [atacs[0],atacs[1],atacs[2]], 
              [adts[0],adts[1],adts[2]]],
    layers = [[rna_input[0],rna_input[1],rna_input[2]], 
              [atac_input[0],atac_input[1],atac_input[2]], 
              [adt_input[0],adt_input[1],adt_input[2]]],
)


del rnas
del atacs
del adts
gc.collect()


new_model = mtg.model.MultiVAE.load_query_data(query,models)
new_model.train(weight_decay=0)
new_model.get_latent_representation()


adata_both = ad.concat([adata, query])
latent_data = adata_both.obsm['latent']

if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:

    latent_df = latent_data


latent_df.to_csv('/home/server/sqy/data/mosaic/result/mapping/multiome/Multigrate.csv', index=False)


# # CITE Transfer

n_batches = 3
rnas = []
adts = []
atacs = []

rna_input = []
adt_input = []
atac_input = []

rna_t = [0,1,2]
adt_t =[0,1,2]
atac_t =[0,1,2]

path = ['CITE/s2d1_cite','CITE/s2d4_cite','CITE/s2d5_cite']


for batch in range(n_batches):
    
    if batch in rna_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/rna.rds"))
        counts_rna = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_rna.shape[0])
        obs['modality'] = np.repeat('CITE', counts_rna.shape[0])
        # np.repeat('rna', counts_rna.shape[0])
        rna_ad = ad.AnnData(counts_rna ,obs = obs)
        rna_ad.obs.index = rds_data[None].keys()
        rna_in = None
        
        # counts_rna = pyreadr.read_r(os.path.join(dir, 'rna' + str(batch + 1) + ".rds")).toarray().T
    else:
        rna_ad = None
        rna_in = None
    
    
    if batch in adt_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/adt.rds"))
        counts_protein = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_protein.shape[0])
        obs['modality'] = np.repeat('CITE', counts_protein.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_protein.shape[0])
        # np.repeat('adt', counts_protein.shape[0])
        adt_ad = ad.AnnData(counts_protein ,obs = obs)
        adt_ad.obs.index = rds_data[None].keys()
        muon.prot.pp.clr(adt_ad)
        adt_in = None
        
    else:
        adt_ad = None
        adt_in = None
    
    if batch in atac_t:
        counts_atac = csr_matrix((counts_rna.shape[0],rds_atac[None].shape[0]))
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_atac.shape[0])
        obs['modality'] = np.repeat('CITE', counts_atac.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_atac.shape[0])
        # np.repeat('atac', counts_atac.shape[0])
        atac_ad = ad.AnnData(counts_atac ,obs = obs)
        atac_ad.obs.index = rna_ad.obs.index
        atac_in = None
        
    else:
        atac_ad = None
        atac_in = None
        
    rnas.append(rna_ad)
    rna_input.append(rna_in)
    
    atacs.append(atac_ad)
    atac_input.append(atac_in)
    
    adts.append(adt_ad)
    adt_input.append(adt_in)


del rds_data
del counts_protein
del counts_rna
import gc
gc.collect()


query = mtg.data.organize_multiome_anndatas(
    adatas = [[rnas[0],rnas[1],rnas[2]], 
              [atacs[0],atacs[1],atacs[2]], 
              [adts[0],adts[1],adts[2]]],
    layers = [[rna_input[0],rna_input[1],rna_input[2]], 
              [atac_input[0],atac_input[1],atac_input[2]], 
              [adt_input[0],adt_input[1],adt_input[2]]],
)

del rnas
del atacs
del adts
gc.collect()


new_model = mtg.model.MultiVAE.load_query_data(query,models)
new_model.train(weight_decay=0)
new_model.get_latent_representation()

adata_both = ad.concat([adata, query])
latent_data = adata_both.obsm['latent']

if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:

    latent_df = latent_data


latent_df.to_csv('/home/server/sqy/data/mosaic/result/mapping/cite/Multigrate.csv', index=False)


# # CITE Multi Transfer


n_batches = 6
rnas = []
atacs = []
adts = []

rna_input = []
atac_input = []
adt_input = []

rna_t = [0,1,2,3,4,5]
atac_t = [3,4,5]
adt_t =[0,1,2]

path = ['CITE/s2d1_cite','CITE/s2d4_cite','CITE/s2d5_cite',
       'Multiome/s2d1_multi','Multiome/s2d4_multi','Multiome/s2d5_multi']


for batch in range(n_batches):
    
    if batch in rna_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/rna.rds"))
        counts_rna = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_rna.shape[0])
        if batch in [0,1,2]:
            obs['modality'] = np.repeat('CITE', counts_rna.shape[0])
        else:
            obs['modality'] = np.repeat('Multiome', counts_rna.shape[0])
        # np.repeat('rna', counts_rna.shape[0])
        rna_ad = ad.AnnData(counts_rna ,obs = obs)
        rna_ad.obs.index = rds_data[None].keys()
        rna_in = None
        
        # counts_rna = pyreadr.read_r(os.path.join(dir, 'rna' + str(batch + 1) + ".rds")).toarray().T
    else:
        rna_ad = None
        rna_in = None
    
    if batch in atac_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/atac.rds"))
        counts_atac = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_atac.shape[0])
        obs['modality'] = np.repeat('Multiome', counts_atac.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_atac.shape[0])
        # np.repeat('atac', counts_atac.shape[0])
        atac_ad = ad.AnnData(counts_atac ,obs = obs)
        atac_ad.obs.index = rds_data[None].keys()
        sc.pp.normalize_total(atac_ad, target_sum=1e4)
        sc.pp.log1p(atac_ad)
        atac_in = None
        
    else:
        atac_ad = None
        atac_in = None
    
    if batch in adt_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/adt.rds"))
        counts_protein = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_protein.shape[0])
        obs['modality'] = np.repeat('CITE', counts_protein.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_protein.shape[0])
        # np.repeat('adt', counts_protein.shape[0])
        adt_ad = ad.AnnData(counts_protein ,obs = obs)
        adt_ad.obs.index = rds_data[None].keys()
        muon.prot.pp.clr(adt_ad)
        adt_in = None
        
    else:
        adt_ad = None
        adt_in = None
        
    rnas.append(rna_ad)
    rna_input.append(rna_in)
    
    atacs.append(atac_ad)
    atac_input.append(atac_in)
    
    adts.append(adt_ad)
    adt_input.append(adt_in)

del rds_data
del counts_protein
del counts_atac
del counts_rna
import gc
gc.collect()


query = mtg.data.organize_multiome_anndatas(
    adatas = [[rnas[0],rnas[1],rnas[2],rnas[3],rnas[4],rnas[5]], 
              [atacs[0],atacs[1],atacs[2],atacs[3],atacs[4],atacs[5]], 
              [adts[0],adts[1],adts[2],adts[3],adts[4],adts[5]]],
    layers = [[rna_input[0],rna_input[1],rna_input[2],rna_input[3],rna_input[4],rna_input[5]], 
              [atac_input[0],atac_input[1],atac_input[2],atac_input[3],atac_input[4],atac_input[5]], 
              [adt_input[0],adt_input[1],adt_input[2],adt_input[3],adt_input[4],adt_input[5]]],
)


del rnas
del atacs
del adts
gc.collect()



new_model = mtg.model.MultiVAE.load_query_data(query,models)
new_model.train(weight_decay=0)
new_model.get_latent_representation()


adata_both = ad.concat([adata, query])
latent_data = adata_both.obsm['latent']

if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:

    latent_df = latent_data


latent_df.to_csv('/home/server/sqy/data/mosaic/result/mapping/cite_multi/Multigrate.csv', index=False)


# # ATAC Transfer


n_batches = 3
atacs = []
rnas = []
adts = []

atac_input = []
rna_input = []
adt_input = []

atac_t = [0,1,2]
rna_t = [0,1,2]
adt_t = [0,1,2]

path = [
       'Multiome/s2d1_multi','Multiome/s2d4_multi','Multiome/s2d5_multi']



for batch in range(n_batches):
    
    if batch in atac_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/atac.rds"))
        counts_atac = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_atac.shape[0])
        obs['modality'] = np.repeat('ATAC', counts_atac.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_atac.shape[0])
        # np.repeat('atac', counts_atac.shape[0])
        atac_ad = ad.AnnData(counts_atac ,obs = obs)
        atac_ad.obs.index = rds_data[None].keys()
        sc.pp.normalize_total(atac_ad, target_sum=1e4)
        sc.pp.log1p(atac_ad)
        atac_in = None
        
    else:
        atac_ad = None
        atac_in = None
    
    if batch in rna_t:
        counts_rna = csr_matrix((counts_atac.shape[0],rds_rna[None].shape[0]))
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_rna.shape[0])
        obs['modality'] = np.repeat('ATAC', counts_rna.shape[0])
        # np.repeat('rna', counts_rna.shape[0])
        rna_ad = ad.AnnData(counts_rna ,obs = obs)
        rna_ad.obs.index = atac_ad.obs.index
        rna_in = None
        
        # counts_rna = pyreadr.read_r(os.path.join(dir, 'rna' + str(batch + 1) + ".rds")).toarray().T
    else:
        rna_ad = None
        rna_in = None
    
    if batch in adt_t:
        counts_protein = csr_matrix((counts_atac.shape[0],rds_adt[None].shape[0]))
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_protein.shape[0])
        obs['modality'] = np.repeat('ATAC', counts_protein.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_protein.shape[0])
        # np.repeat('adt', counts_protein.shape[0])
        adt_ad = ad.AnnData(counts_protein ,obs = obs)
        adt_ad.obs.index = atac_ad.obs.index
        adt_in = None
        
    else:
        adt_ad = None
        adt_in = None
        
    rnas.append(rna_ad)
    rna_input.append(rna_in)
    
    atacs.append(atac_ad)
    atac_input.append(atac_in)
    
    adts.append(adt_ad)
    adt_input.append(adt_in)


del rds_data
del counts_protein
del counts_atac
del counts_rna
import gc
gc.collect()


query = mtg.data.organize_multiome_anndatas(
    adatas = [[rnas[0],rnas[1],rnas[2]], 
              [atacs[0],atacs[1],atacs[2]], 
              [adts[0],adts[1],adts[2]]],
    layers = [[rna_input[0],rna_input[1],rna_input[2]], 
              [atac_input[0],atac_input[1],atac_input[2]], 
              [adt_input[0],adt_input[1],adt_input[2]]],
)

del rnas
del atacs
del adts
gc.collect()



new_model = mtg.model.MultiVAE.load_query_data(query,models)
new_model.train(weight_decay=0)
new_model.get_latent_representation()


adata_both = ad.concat([adata, query])
latent_data = adata_both.obsm['latent']

if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:

    latent_df = latent_data


latent_df.to_csv('/home/server/sqy/data/mosaic/result/mapping/atac/Multigrate.csv', index=False)


# # ADT ATAC Transfer

n_batches = 6
atacs = []
adts = []
rnas = []

atac_input = []
adt_input = []
rna_input = []

atac_t = [3,4,5]
adt_t =[0,1,2]
rna_t =[0,1,2]

path = ['CITE/s2d1_cite','CITE/s2d4_cite','CITE/s2d5_cite',
       'Multiome/s2d1_multi','Multiome/s2d4_multi','Multiome/s2d5_multi']


for batch in range(n_batches):
    
    if batch in atac_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/atac.rds"))
        counts_atac = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_atac.shape[0])
        obs['modality'] = np.repeat('ATAC', counts_atac.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_atac.shape[0])
        # np.repeat('atac', counts_atac.shape[0])
        atac_ad = ad.AnnData(counts_atac ,obs = obs)
        atac_ad.obs.index = rds_data[None].keys()
        sc.pp.normalize_total(atac_ad, target_sum=1e4)
        sc.pp.log1p(atac_ad)
        atac_in = None
        
    else:
        atac_ad = None
        atac_in = None
    
    if batch in adt_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/adt.rds"))
        counts_protein = csr_matrix(np.array(rds_data[None]).T)
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_protein.shape[0])
        obs['modality'] = np.repeat('ADT', counts_protein.shape[0])
        # obs['modality'] = np.repeat('batch' + str(batch + 1), counts_protein.shape[0])
        # np.repeat('adt', counts_protein.shape[0])
        adt_ad = ad.AnnData(counts_protein ,obs = obs)
        adt_ad.obs.index = rds_data[None].keys()
        muon.prot.pp.clr(adt_ad)
        adt_in = None
        
    else:
        adt_ad = None
        adt_in = None
    
    if batch in rna_t:
        counts_rna = csr_matrix((counts_protein.shape[0],rds_rna[None].shape[0]))
        obs = pd.DataFrame()
        obs["samplename"] = np.repeat('batch' + str(batch+6), counts_rna.shape[0])
        obs['modality'] = np.repeat('ATAC', counts_rna.shape[0])
        # np.repeat('rna', counts_rna.shape[0])
        rna_ad = ad.AnnData(counts_rna ,obs = obs)
        rna_ad.obs.index = adt_ad.obs.index
        rna_in = None
        
        # counts_rna = pyreadr.read_r(os.path.join(dir, 'rna' + str(batch + 1) + ".rds")).toarray().T
    else:
        rna_ad = None
        rna_in = None
        
    rnas.append(rna_ad)
    rna_input.append(rna_in)
    
    atacs.append(atac_ad)
    atac_input.append(atac_in)
    
    adts.append(adt_ad)
    adt_input.append(adt_in)


del rds_data
del counts_protein
del counts_atac
del counts_rna
import gc
gc.collect()



query = mtg.data.organize_multiome_anndatas(
    adatas = [[rnas[0],rnas[1],rnas[2],rnas[3],rnas[4],rnas[5]], 
              [atacs[0],atacs[1],atacs[2],atacs[3],atacs[4],atacs[5]], 
              [adts[0],adts[1],adts[2],adts[3],adts[4],adts[5]]],
    layers = [[rna_input[0],rna_input[1],rna_input[2],rna_input[3],rna_input[4],rna_input[5]], 
              [atac_input[0],atac_input[1],atac_input[2],atac_input[3],atac_input[4],atac_input[5]], 
              [adt_input[0],adt_input[1],adt_input[2],adt_input[3],adt_input[4],adt_input[5]]],
)


del rnas
del atacs
del adts
gc.collect()


new_model = mtg.model.MultiVAE.load_query_data(query,models)
new_model.train(weight_decay=0)
new_model.get_latent_representation()


adata_both = ad.concat([adata, query])
latent_data = adata_both.obsm['latent']

if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:

    latent_df = latent_data


latent_df.to_csv('/home/server/sqy/data/mosaic/result/mapping/adt_atac/Multigrate.csv', index=False)

