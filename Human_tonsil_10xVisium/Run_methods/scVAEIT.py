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

dir = "./10x_Visium_human_tonsil/"
n_batches = 3

batch_name = None

if batch_name != None and len(batch_name) != n_batches :
    raise ValueError("batch_name wrong !")

rnas = []
atacs = []
adts = []

rna_t = [0,1]
adt_t =[0,2]

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
path = ['slice1','slice2','slice3']

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
    
    if batch in adt_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/adt.rds"))
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
    adts.append(counts_protein)
    
# batch_idx = [item for vector in batch_id for item in vector]
batches = np.concatenate(batch_id, axis=0)
new_rna = []
new_adt = []

for i in range(n_batches):
    if np.all(rnas[i] == None):
        new_rna.append(np.zeros((cell_num[i],gene_num)))
        # new_rna[i] = np.zeros(cell_num[i],gene_num)
    else:
        new_rna.append(rnas[i])
    
    if np.all(adts[i] == None):
        new_adt.append(np.zeros((cell_num[i],protein_num)))
        #new_adt[i] = np.zeros(cell_num[i],protein_num)
    else:
        new_adt.append(adts[i])

rna_mat = np.concatenate(new_rna, axis=0,dtype = np.float32)
adt_mat = np.concatenate(new_adt, axis=0,dtype = np.float32)

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

dim_input_arr = np.array([gene_num, protein_num])
masks = - np.ones((n_batches, np.sum(dim_input_arr)), dtype=np.float32)
for i in range(n_batches):
    if np.any(rnas[i] != None):
        tmp = np.arange(0,gene_num,1)
        masks[i,tmp] = 0.
    
    if np.any(adts[i] != None):
        tmp = np.arange(gene_num,gene_num+protein_num,1)
        masks[i,tmp] = 0. 

masks = tf.convert_to_tensor(masks, dtype=tf.float32)

for i in rna_t:
    rna_mat1[batches[:,-1]==i,:] = np.log(rna_mat1[batches[:,-1]==i,:]/np.sum(rna_mat1[batches[:,-1]==i,:], axis=1, keepdims=True)*1e4+1.)

for j in adt_t:
    adt_mat1[batches[:,-1]==j,:] = np.log(adt_mat1[batches[:,-1]==j,:]/np.sum(adt_mat1[batches[:,-1]==j,:], axis=1, keepdims=True)*1e4+1.)
# adt_mat = np.log(adt_mat/np.sum(adt_mat, axis=1, keepdims=True)*1e4+1.)

data = np.c_[rna_mat1, adt_mat1]

import sys # add parent folder path if not installed via PyPI
import scVAEIT
print(scVAEIT.__version__)

from scVAEIT.VAEIT import VAEIT
path_root = '/home/server/sqy/code/scVAEIT/result/'

config = {
    # Dimension of input features for [RNA, ADT, peaks]
    'dim_input_arr': dim_input_arr,

    # Blocks for [RNA, ADT, peaks at chr1, ... peaks at chr22]
    'dim_block': [gene_num, protein_num], # input dimension of blocks
    'dist_block':['NB','NB'] , # distributions of blocks
    'dim_block_enc':np.array([256, 128]), # dimension of first layer of the encoder
    'dim_block_dec':np.array([256, 128]), # dimension of first layer of the decoder
    'block_names':np.array(['rna', 'adt']), # names of blocks
    'uni_block_names':np.array(['rna','adt']), # names for modalities
    'dim_block_embed':np.array([16, 8])*2, # mask embedding dimension

    # Internal network structure
    'dimensions':[256], # dimension of latent layers of encoder; the reversed is used for decoder
    'dim_latent':32, # the latent dimension bewteen the encoder and decoder
    
    # Weights
    'beta_unobs':2./3., # weight for masked out observation; weight for observerd values will be 1-beta_unobs.
    'beta_modal':np.array([0.15,0.85]), # weights for 3 modalities

    # Masking probability
    "p_feat" : 0.2, # probablity of randomly masking out an entry
    "p_modal" : np.ones(2)/2, # probability of masking out the 3 modalities
    
    'gamma':lambda epoch: 2 * 0.8 ** (epoch / 50) # MMD loss to correct strong batch effects of multiple datasets
}
del rnas
del adts
del rna_t
del adt_t
del gene_id
del protein_id
del rds_data
del counts_rna
del counts_protein
del new_rna
del new_adt
del rna_mat
del rna_mat1
del adt_mat
del adt_mat1
import gc
gc.collect()

# batches[:,0] = 0.
batch = batches[:,0]
# batch = batch.reshape((batch.shape[0],1))

model = VAEIT(config, data,
              masks,batch,
              batches,None,None)
hist = model.train(
        num_epoch=300, batch_size=100, save_every_epoch=300,
        verbose=True, checkpoint_dir=path_root+'checkpoint/')

latent =  model.get_latent_z()
# np.savetxt('/home/server/sqy/data/mosaic/result/bench/ABseq/scVAEIT.csv', latent, delimiter=',')

if isinstance(latent, np.ndarray):
    latent = latent[:, :20]
    latent_df = pd.DataFrame(latent)
else:
    latent_df = latent.iloc[:, :20]

latent_df.to_csv('./10xVisium_tonsil/scVAEIT.csv', index=False)

