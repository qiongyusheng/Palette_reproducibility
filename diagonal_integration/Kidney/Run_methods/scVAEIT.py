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

dir = "./"
n_batches = 2

batch_name = None

if batch_name != None and len(batch_name) != n_batches :
    raise ValueError("batch_name wrong !")

rnas = []
atacs = []

rna_t = [0]
atac_t = [1]

#用于判断是否已经保存了特征名称
gene_id = False
peak_id = False

cell_num = [None]*n_batches
gene_num = None
peak_num = None

#batch_id = [[] for _ in range(n_batches)]
batch_id = [None]*n_batches
path = ['RNA','ATAC']

for batch in range(n_batches):
    
    if batch in rna_t:
        rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/X.rds"))
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
        if peak_id == False:
            peak_name = np.array(rds_data[None].index, dtype=str)#rds_data[None].index
            peak_id == True

        # rds_data = pyreadr.read_r(os.path.join(dir, path[batch] + "/atac.rds"))
        # counts_atac = np.array(rds_data[None]).T
        # counts_atac = counts_atac#[:, 1:80000]
        
        # if peak_id == False:
        #     peak_name = np.array(rds_data[None].index, dtype=str)#[1:80000] #rds_data[None].index
        #     peak_id = True
            
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
peak_name1 = np.array(peak_name[non_zero_col_indices],dtype=str)
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

# 1. 首先计算最终矩阵的形状
final_cols = rna_mat1.shape[1] + atac_mat1.shape[1]
rows = rna_mat1.shape[0]

# 2. 预先分配内存（使用float32来节省内存）
data = np.empty((rows, final_cols), dtype=np.float32)

# 3. 分段填充数据
print("Starting to fill data matrix...")

# 填充rna_mat1
print("Copying rna_mat1...")
col_start = 0
col_end = rna_mat1.shape[1]
data[:, col_start:col_end] = rna_mat1
del rna_mat1  # 立即释放内存
print("rna_mat1 copied")

# 分批填充atac_mat1（因为它最大）
print("Copying atac_mat1...")
col_start = col_end
col_end = col_start + atac_mat1.shape[1]

# 使用更小的批次
batch_size = 1000  # 每次处理1000列
for i in range(0, atac_mat1.shape[1], batch_size):
    end_idx = min(i + batch_size, atac_mat1.shape[1])
    data[:, col_start+i:col_start+end_idx] = atac_mat1[:, i:end_idx]
    print(f"Processed columns {i} to {end_idx} of atac_mat1")

del atac_mat1  # 释放内存
print("atac_mat1 copied")

print("Data matrix creation completed!")

import sys # add parent folder path if not installed via PyPI
import scVAEIT
print(scVAEIT.__version__)

from scVAEIT.VAEIT import VAEIT
path_root = '/home/server/sqy/code/scVAEIT/result/'

config = {
    # Dimension of input features for [RNA, ADT, peaks]
    'dim_input_arr': dim_input_arr,

    # Blocks for [RNA, ADT, peaks at chr1, ... peaks at chr22]
    'dim_block': np.append([gene_num], chunk_atac), # input dimension of blocks
    'dist_block':['NB'] + ['Bernoulli' for _ in chunk_atac], # distributions of blocks
    'dim_block_enc':np.array([256] + [16 for _ in chunk_atac]), # dimension of first layer of the encoder
    'dim_block_dec':np.array([256] + [16 for _ in chunk_atac]), # dimension of first layer of the decoder
    'block_names':np.array(['rna'] + ['atac' for _ in range(len(chunk_atac))]), # names of blocks
    'uni_block_names':np.array(['rna','atac']), # names for modalities
    'dim_block_embed':np.array([16] + [1 for _ in range(len(chunk_atac))])*2, # mask embedding dimension

    # Internal network structure
    'dimensions':[256], # dimension of latent layers of encoder; the reversed is used for decoder
    'dim_latent':32, # the latent dimension bewteen the encoder and decoder
    
    # Weights
    'beta_unobs':2./3., # weight for masked out observation; weight for observerd values will be 1-beta_unobs.
    'beta_modal':np.array([0.93,0.07]), # weights for 3 modalities

    # Masking probability
    "p_feat" : 0.2, # probablity of randomly masking out an entry
    "p_modal" : np.ones(2)/2, # probability of masking out the 3 modalities
    
    'gamma':lambda epoch: 2 * 0.8 ** (epoch / 50) # MMD loss to correct strong batch effects of multiple datasets
}

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

if isinstance(latent, np.ndarray):
    latent = latent
    latent_df = pd.DataFrame(latent)
else:
    latent_df = latent.iloc

latent_df.to_csv('./kidney/scVAEIT.csv', index=False)
