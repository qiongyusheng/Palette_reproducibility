import os
import sys
os.environ["OMP_NUM_THREADS"] = "11"
os.environ["OPENBLAS_NUM_THREADS"] = "8" # export OPENBLAS_NUM_THREADS=4 
os.environ["MKL_NUM_THREADS"] = "11" # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = "8" # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = "11" # export NUMEXPR_NUM_THREADS=6
# os.environ["NUMBA_CACHE_DIR"]='/tmp/numba_cache'
import numpy as np
import pandas as pd
import scipy as sp
import scipy.sparse
import h5py
import pyreadr
# import anndata as ad#不必要
# import muon

import tensorflow as tf
import matplotlib.pyplot as plt
import scanpy as sc

import os, sys, time
import numpy as np
import pandas as pd
import scipy as sp
import scipy.sparse
import h5py
import tensorflow as tf
import tensorflow_probability as tfp
import scanpy as sc
from tensorflow.keras.layers import Layer, Dense, BatchNormalization, LeakyReLU
from tensorflow.keras.utils import Progbar
tfd = tfp.distributions

import matplotlib.pyplot as plt
from types import SimpleNamespace
from sklearn.model_selection import train_test_split

from types import SimpleNamespace

from scipy.sparse import csr_matrix, coo_matrix, csc_matrix,lil_matrix
from scipy.io import mmread, mmwrite
sys.path.insert(0, '/home/server/sqy/code/') 
from scVAEIT.VAEIT import scVAEIT
path_root = '/home/server/sqy/code/scVAEIT/result/'

dir = "./PBMC_TEA/"
n_batches = 4

batch_name = None

if batch_name != None and len(batch_name) != n_batches :
    raise ValueError("batch_name wrong !")

rnas = []
atacs = []
adts = []

rna_t = [1,2]
atac_t = [0,1]
adt_t =[2,3]

#用于判断是否已经保存了特征名称
gene_id = False
peak_id = False
protein_id = False

cell_num = [None]*n_batches
gene_num = None
peak_num = None
protein_num = None

#batch_id = [[] for _ in range(n_batches)]
batch_id = [None]*n_batches

for batch in range(n_batches):
    
    if batch in rna_t:
        rds_data = pyreadr.read_r(os.path.join(dir, 'B' + str(batch + 2) + "/rna.rds"))
        counts_rna = np.array(rds_data[None]).T
        if gene_id == False:
            gene_name = np.array(rds_data[None].index, dtype=str)#rds_data[None].index
            gene_id == True
            
        if cell_num[batch] == None:
            cell_num[batch] = counts_rna.shape[0]
        
        if gene_num == None:
            gene_num = counts_rna.shape[1]
            
        # if len(batch_id[batch]) == 0:
        if np.all(batch_id[batch] == None):
            
            batch_id[batch] = np.full((counts_rna.shape[0],2), batch, dtype = np.float32)
            
            #if batch_name == None:
            #    batch_id[batch] = ['batch' + str(batch + 1)] * counts_rna.shape[0]
            #else:
            #    batch_id[batch] = batch_name[batch] * counts_rna.shape[0]
        
    else:
        counts_rna = None
    
    if batch in atac_t:
        rds_data = pyreadr.read_r(os.path.join(dir, 'B' + str(batch + 2) + "/atac.rds"))
        counts_atac = np.array(rds_data[None]).T
        counts_atac = counts_atac#[:, 1:80000]
        
        if peak_id == False:
            peak_name = np.array(rds_data[None].index, dtype=str)#[1:80000] #rds_data[None].index
            peak_id = True
            
        if cell_num[batch] == None:
            cell_num[batch] = counts_atac.shape[0]
        
        if peak_num == None:
            peak_num = counts_atac.shape[1] 
        
        #if len(batch_id[batch]) == 0:
        if np.all(batch_id[batch] == None):
            
            batch_id[batch] = np.full((counts_atac.shape[0],2), batch, dtype = np.float32)
            
            #if batch_name == None:
            #    batch_id[batch] = ['batch' + str(batch + 1)] * counts_atac.shape[0]
            #else:
            #    batch_id[batch] = batch_name[batch] * counts_atac.shape[0]
        
        
    else:
        counts_atac = None
    
    if batch in adt_t:
        rds_data = pyreadr.read_r(os.path.join(dir, 'B' + str(batch + 2) + "/adt.rds"))
        counts_protein = np.array(rds_data[None]).T
        
        if protein_id == False:
            protein_name = np.array(rds_data[None].index, dtype=str) # rds_data[None].index
            protein_id = True
            
        if cell_num[batch] == None:
            cell_num[batch] = counts_protein.shape[0]
        
        if protein_num == None:
            protein_num = counts_protein.shape[1]
        
        #if len(batch_id[batch]) == 0:
        if np.all(batch_id[batch] == None):
            
            batch_id[batch] = np.full((counts_protein.shape[0],2), batch, dtype = np.float32)
            
            #if batch_name == None:
            #    batch_id[batch] = ['batch' + str(batch + 1)] * counts_protein.shape[0]
            #else:
            #    batch_id[batch] = batch_name[batch] * counts_protein.shape[0]
        
    else:
        counts_protein = None
        
        
    rnas.append(counts_rna)
    atacs.append(counts_atac)
    adts.append(counts_protein)
    
# batch_idx = [item for vector in batch_id for item in vector]
batches = np.concatenate(batch_id, axis=0)

import scipy.sparse as sp
tmp = sp.vstack([atacs[0], atacs[3]])
import anndata as ad
import scanpy as sc
adata = ad.AnnData(tmp,dtype=np.float32)
adata.obs['batch'] = ["batch1"] * atacs[0].shape[0] + ["batch2"] * atacs[1].shape[0]
# adata.layers['counts'] = adata.X.copy()
adata.var_names = peak_name
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=50000, batch_key='batch')

adata_hvg = adata[:, adata.var.highly_variable].copy()
regions_new = adata_hvg.var_names

adata_hvg0 = adata_hvg[adata.obs['batch'] == 'batch1',:]
hvg0_d = np.array(adata_hvg0.X.todense())
atacs[0] = hvg0_d
del adata_hvg0
del hvg0_d

adata_hvg1 = adata_hvg[adata.obs['batch'] == 'batch2',:]
hvg1_d = np.array(adata_hvg1.X.todense())
# counts_atac1 = utils.preprocess(hvg1_d, modality = "ATAC")
atacs[1] = hvg1_d
del adata_hvg1
# del counts_atac1
del hvg1_d

del tmp
del adata
del adata_hvg

peak_num = atacs[0].shape[1]
new_rna = []
new_atac = []
new_adt = []

for i in range(n_batches):
    if np.all(rnas[i] == None):
        new_rna.append(np.zeros((cell_num[i],gene_num)))
        # new_rna[i] = np.zeros(cell_num[i],gene_num)
    else:
        new_rna.append(rnas[i])
    
    if np.all(atacs[i] == None):
        new_atac.append(np.zeros((cell_num[i],peak_num)))
        #new_atac[i] = np.zeros(cell_num[i],peak_num)
    else:
        new_atac.append(atacs[i])
    
    if np.all(adts[i] == None):
        new_adt.append(np.zeros((cell_num[i],protein_num)))
        #new_adt[i] = np.zeros(cell_num[i],protein_num)
    else:
        new_adt.append(adts[i])

rna_mat = np.concatenate(new_rna, axis=0,dtype = np.float32)
atac_mat = np.concatenate(new_atac, axis=0,dtype = np.float32)
adt_mat = np.concatenate(new_adt, axis=0,dtype = np.float32)

col_sums = np.sum(atac_mat, axis=0)
non_zero_col_indices = col_sums != 0
atac_mat1 = atac_mat[:,non_zero_col_indices]
peak_name1 = peak_name[non_zero_col_indices]
peak_num = atac_mat1.shape[1]
col_sums = np.sum(rna_mat, axis=0)
non_zero_col_indices = col_sums != 0
rna_mat1 = rna_mat[:,non_zero_col_indices]
gene_name1 = gene_name[non_zero_col_indices]
gene_num = rna_mat1.shape[1]
col_sums = np.sum(adt_mat, axis=0)
non_zero_col_indices = col_sums != 0
adt_mat1 = adt_mat[:,non_zero_col_indices]
protein_name1 =  protein_name[non_zero_col_indices]
protein_num = adt_mat1.shape[1]

dim_input_arr = np.array([gene_num, protein_num, peak_num])
masks = - np.ones((n_batches, np.sum(dim_input_arr)), dtype=np.float32)
for i in range(n_batches):
    if np.any(rnas[i] != None):
        tmp = np.arange(0,gene_num,1)
        masks[i,tmp] = 0.
    
    if np.any(adts[i] != None):
        tmp = np.arange(gene_num,gene_num+protein_num,1)
        masks[i,tmp] = 0. 
    
    if np.any(atacs[i] != None):
        tmp = np.arange(gene_num+protein_num, gene_num+protein_num+peak_num, 1)
        masks[i,tmp] = 0.

masks = tf.convert_to_tensor(masks, dtype=tf.float32)
# Count peaks in each chromosome (assuming they are ordered)
chunk_atac = np.array([
    np.sum(np.char.startswith(peak_name1, 'chr%d-'%i)) for i in range(1,23)], dtype=np.int32)
for i in rna_t:
    rna_mat1[batches[:,-1]==i,:] = np.log(rna_mat1[batches[:,-1]==i,:]/np.sum(rna_mat1[batches[:,-1]==i,:], axis=1, keepdims=True)*1e4+1.)
#rna_mat = np.log(rna_mat/np.sum(rna_mat, axis=1, keepdims=True)*1e4+1.)
for j in adt_t:
    adt_mat1[batches[:,-1]==j,:] = np.log(adt_mat1[batches[:,-1]==j,:]/np.sum(adt_mat1[batches[:,-1]==j,:], axis=1, keepdims=True)*1e4+1.)
# adt_mat = np.log(adt_mat/np.sum(adt_mat, axis=1, keepdims=True)*1e4+1.)
atac_mat1[atac_mat1>0.] = 1.

data = np.c_[rna_mat1, adt_mat1, atac_mat1]

config = {
    'dim_input_arr': dim_input_arr,
    'dimensions':[256], 
    'dim_latent':32,
    'dim_block': np.append([gene_num, protein_num], chunk_atac), # input dimension of blocks
    'dist_block':['NB','NB'] + ['Bernoulli' for _ in chunk_atac], 
    'dim_block_enc':np.array([256, 128] + [16 for _ in chunk_atac]),
    'dim_block_dec':np.array([256, 128] + [16 for _ in chunk_atac]),
    'block_names':np.array(['rna', 'adt'] + ['atac' for _ in range(len(chunk_atac))]),
    'uni_block_names':np.array(['rna','adt','atac']),
    'dim_block_embed':np.array([16, 8] + [1 for _ in range(len(chunk_atac))]),

    'beta_kl':2.,
    'beta_unobs':2./3.,
    'beta_modal':np.array([0.14,0.85,0.01]),
    'beta_reverse':0.2,

    "p_feat" : 0.2,
    "p_modal" : np.ones(3)/3,
    
}
config = SimpleNamespace(**config)
del rnas
del atacs
del adts
del rna_t
del atac_t
del adt_t
del gene_id
del peak_id
del protein_id
del rds_data
del counts_rna
del counts_atac
del counts_protein
del new_rna
del new_atac
del new_adt
del rna_mat
del rna_mat1
del atac_mat
del atac_mat1
del adt_mat
del adt_mat1
import gc
gc.collect()
# batches[:,0] = 0.
batch = batches[:,0]
batch = batch.reshape((batch.shape[0],1))

model = scVAEIT(config, data,
              masks, 
              batch)

hist = model.train(
        valid = False, num_epoch = 300, batch_size = 521, save_every_epoch = 300,
        verbose = True, checkpoint_dir = None)
latent =  model.get_latent_z()

if isinstance(latent, np.ndarray):
    latent = latent[:, :20]
    latent_df = pd.DataFrame(latent)
else:
    latent_df = latent.iloc[:, :20]

latent_df.to_csv('./Benchmarking/Unsupervised/TEA_s2/output/scVAEIT.csv', index=False)