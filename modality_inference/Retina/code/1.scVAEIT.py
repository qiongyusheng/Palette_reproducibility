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

dir = "./imputation/input/Retina/"
n_batches = 8

batch_name = None

if batch_name != None and len(batch_name) != n_batches :
    raise ValueError("batch_name wrong !")

rnas = []
atacs = []

rna_t = [0,1,2,3,4]
atac_t = [0,1,5,6,7]

#用于判断是否已经保存了特征名称
gene_id = False
peak_id = False

cell_num = [None]*n_batches
gene_num = None
peak_num = None

#batch_id = [[] for _ in range(n_batches)]
batch_id = [None]*n_batches
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
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/atac.rds"))
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
      
    rnas.append(counts_rna)
    atacs.append(counts_atac)


batches = np.concatenate(batch_id, axis=0)

batches = np.concatenate(batch_id, axis=0)

new_rna = []
new_atac = []

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

rna_mat = np.concatenate(new_rna, axis=0,dtype = np.float32)
atac_mat = np.concatenate(new_atac, axis=0,dtype = np.float32)

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

dim_input_arr = np.array([gene_num, peak_num])
masks = - np.ones((n_batches, np.sum(dim_input_arr)), dtype=np.float32)

for i in range(n_batches):
    if np.any(rnas[i] != None):
        tmp = np.arange(0,gene_num,1)
        masks[i,tmp] = 0.
    
    if np.any(atacs[i] != None):
        tmp = np.arange(gene_num ,gene_num+peak_num, 1)
        masks[i,tmp] = 0.

masks = tf.convert_to_tensor(masks, dtype=tf.float32)

chunk_atac = np.array([
    np.sum(np.char.startswith(peak_name1, 'chr%d-'%i)) for i in range(1,23)], dtype=np.int32)

for i in rna_t:
    rna_mat1[batches[:,-1]==i,:] = np.log(rna_mat1[batches[:,-1]==i,:]/np.sum(rna_mat1[batches[:,-1]==i,:], axis=1, keepdims=True)*1e4+1.)
atac_mat1[atac_mat1>0.] = 1.

data = np.c_[rna_mat1, atac_mat1]

config = {
    'dim_input_arr': dim_input_arr,
    'dimensions':[256], 
    'dim_latent':32,
    'dim_block': np.append([gene_num], chunk_atac), # input dimension of blocks
    'dist_block':['NB'] + ['Bernoulli' for _ in chunk_atac], 
    'dim_block_enc':np.array([256] + [16 for _ in chunk_atac]),
    'dim_block_dec':np.array([256] + [16 for _ in chunk_atac]),
    'block_names':np.array(['rna'] + ['atac' for _ in range(len(chunk_atac))]),
    'uni_block_names':np.array(['rna','atac']),
    'dim_block_embed':np.array([32] + [2 for _ in range(len(chunk_atac))]),

    'beta_kl':2.,
    'beta_unobs':2./3.,
    'beta_modal':np.array([0.93,0.07]),
    'beta_reverse':0.2,

    "p_feat" : 0.2,
    "p_modal" : np.ones(2)/2,
    
}
config = SimpleNamespace(**config)

del rnas
del atacs
del rna_t
del atac_t
del gene_id
del peak_id
del rds_data
del counts_rna
del counts_atac
del new_rna
del new_atac
del rna_mat
del rna_mat1
del atac_mat
del atac_mat1
import gc
gc.collect()

batch = batches[:,0]
batch = batch.reshape((batch.shape[0],1))

model = scVAEIT(config, data,
              masks, 
              batch)

hist = model.train(
        valid = False, num_epoch = 300, batch_size = 100, save_every_epoch = 300,
        verbose = True, checkpoint_dir = None)

denoised_data = model.get_denoised_data(
    batch_size_inference=128, L=50)
batch = batches[:,0]

# LGS2OD atac
LGS2OD_atac = denoised_data[batch == 2.,model.vae.config.dim_input_arr[0]:(model.vae.config.dim_input_arr[0]+model.vae.config.dim_input_arr[1])]


if isinstance(LGS2OD_atac, np.ndarray):
    latent_df = pd.DataFrame(LGS2OD_atac)
else:
    
    latent_df = LGS2OD_atac


latent_df.to_csv('./imputation/Retina/scVAEIT/LGS2OD_atac.csv', index=False)

del LGS2OD_atac
del latent_df

# LGS2OS atac
LGS2OS_atac = denoised_data[batch == 3.,model.vae.config.dim_input_arr[0]:(model.vae.config.dim_input_arr[0]+model.vae.config.dim_input_arr[1])]


if isinstance(LGS2OS_atac, np.ndarray):
    latent_df = pd.DataFrame(LGS2OS_atac)
else:
    
    latent_df = LGS2OS_atac


latent_df.to_csv('./imputation/Retina/scVAEIT/LGS2OS_atac.csv', index=False)

del LGS2OS_atac
del latent_df

# LGS3OD atac
LGS3OD_atac = denoised_data[batch == 4.,model.vae.config.dim_input_arr[0]:(model.vae.config.dim_input_arr[0]+model.vae.config.dim_input_arr[1])]


if isinstance(LGS3OD_atac, np.ndarray):
    latent_df = pd.DataFrame(LGS3OD_atac)
else:
    
    latent_df = LGS3OD_atac


latent_df.to_csv('./imputation/Retina/scVAEIT/LGS3OD_atac.csv', index=False)

del LGS3OD_atac
del latent_df

# LGS3OS rna
LGS3OS_rna = denoised_data[batch == 5.,0:model.vae.config.dim_input_arr[0]]


if isinstance(LGS3OS_rna, np.ndarray):
    latent_df = pd.DataFrame(LGS3OS_rna)
else:
    
    latent_df = LGS3OS_rna


latent_df.to_csv('./imputation/Retina/scVAEIT/LGS3OS_rna.csv', index=False)

del LGS3OS_rna
del latent_df

# LVG1OD rna
LVG1OD_rna = denoised_data[batch == 6.,0:model.vae.config.dim_input_arr[0]]


if isinstance(LVG1OD_rna, np.ndarray):
    latent_df = pd.DataFrame(LVG1OD_rna)
else:
    
    latent_df = LVG1OD_rna


latent_df.to_csv('./imputation/Retina/scVAEIT/LVG1OD_rna.csv', index=False)

del LVG1OD_rna
del latent_df

# LVG1OS rna
LVG1OS_rna = denoised_data[batch == 7.,0:model.vae.config.dim_input_arr[0]]


if isinstance(LVG1OS_rna, np.ndarray):
    latent_df = pd.DataFrame(LVG1OS_rna)
else:
    
    latent_df = LVG1OS_rna


latent_df.to_csv('./imputation/Retina/scVAEIT/LVG1OS_rna.csv', index=False)

del LVG1OS_rna
del latent_df

peak = pd.DataFrame(peak_name1)
peak.to_csv('./imputation/Retina/scVAEIT/peak.csv', index=False)

gene = pd.DataFrame(gene_name1)
gene.to_csv('./imputation/Retina/scVAEIT/gene.csv', index=False)