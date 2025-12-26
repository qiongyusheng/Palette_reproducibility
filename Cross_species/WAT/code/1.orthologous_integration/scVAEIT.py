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

dir = "./human_mouse_WAT/"
n_batches = 36

batch_name = None

if batch_name != None and len(batch_name) != n_batches :
    raise ValueError("batch_name wrong !")

rnas = []

rna_t = range(n_batches)
gene_id = False

cell_num = [None]*n_batches
gene_num = None
batch_id = [None]*n_batches
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
        
        
    rnas.append(counts_rna)
    
# batch_idx = [item for vector in batch_id for item in vector]
batches = np.concatenate(batch_id, axis=0)

new_rna = []

for i in range(n_batches):
    if np.all(rnas[i] == None):
        new_rna.append(np.zeros((cell_num[i],gene_num)))
        # new_rna[i] = np.zeros(cell_num[i],gene_num)
    else:
        new_rna.append(rnas[i])

rna_mat = np.concatenate(new_rna, axis=0,dtype = np.float32)

col_sums = np.sum(rna_mat, axis=0)
non_zero_col_indices = col_sums != 0
rna_mat1 = rna_mat[:,non_zero_col_indices]
gene_name1 = gene_name[non_zero_col_indices]
gene_num = rna_mat1.shape[1]

dim_input_arr = np.array([gene_num])
masks = - np.ones((n_batches, np.sum(dim_input_arr)), dtype=np.float32)
for i in range(n_batches):
    if np.any(rnas[i] != None):
        tmp = np.arange(0,gene_num,1)
        masks[i,tmp] = 0.

masks = tf.convert_to_tensor(masks, dtype=tf.float32)

for i in rna_t:
    rna_mat1[batches[:,-1]==i,:] = np.log(rna_mat1[batches[:,-1]==i,:]/np.sum(rna_mat1[batches[:,-1]==i,:], axis=1, keepdims=True)*1e4+1.)

data = rna_mat1


config = {
    'dim_input_arr': dim_input_arr,
    'dimensions':[256], 
    'dim_latent':32,
    'dim_block': [gene_num], # input dimension of blocks
    'dist_block':['NB'] , 
    'dim_block_enc':np.array([256]),
    'dim_block_dec':np.array([256]),
    'block_names':np.array(['rna']),
    'uni_block_names':np.array(['rna']),
    'dim_block_embed':np.array([16]),

    'beta_kl':2.,
    'beta_unobs':2./3.,
    'beta_modal':np.array([1]),
    'beta_reverse':0.2,

    "p_feat" : 0.2,
    "p_modal" : np.ones(1)/1,
    
}
config = SimpleNamespace(**config)

del rnas
del rna_t
del gene_id
del rds_data
del counts_rna
del new_rna
del rna_mat
del rna_mat1
import gc
gc.collect()


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

latent_df.to_csv('./scVAEIT_homo.csv', index=False)