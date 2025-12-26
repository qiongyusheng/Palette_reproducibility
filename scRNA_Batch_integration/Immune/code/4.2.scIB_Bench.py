import scanpy as sc
import numpy as np
import pandas as pd
import pyreadr
import anndata
import warnings
warnings.filterwarnings('ignore')
from scib_metrics.benchmark import Benchmarker, BatchCorrection, BioConservation


# In[2]:


out_dir = './'
meta_dir = './'

method_mapping = {'Palette Unsupervised':"Palette_unsupervised",
                  'Palette Supervised':"Palette_supervised",
                  'STACAS':'STACAS',
                  'ssSTACAS':'ssSTACAS',
                  'SIGNAL':'SIGNAL',
                  'scPoli':'scpoli_',
                  'scANVI':'scanvi_',
                  'fastMNN':'fastMNN',
                  'Harmony':'Harmony',
                  'Seurat':'Seurat',
                  'RAW':'raw',
                  
                 'Palette miss 5%':'label_miss/Palette_miss5',
                 'Palette miss 10%':'label_miss/Palette_miss10',
                 'Palette miss %':'label_miss/Palette_miss',
                 'SIGNAL miss 5%':'label_miss/SIGNAL_miss5',
                 'SIGNAL miss 10%':'label_miss/SIGNAL_miss10',
                 'SIGNAL miss %':'label_miss/SIGNAL_miss',
                 'ssSTACAS miss 5%':'label_miss/ssSTACAS_miss5',
                 'ssSTACAS miss 10%':'label_miss/ssSTACAS_miss10',
                 'ssSTACAS miss %':'label_miss/ssSTACAS_miss',
                 'scANVI miss 5%':'label_miss/scanvi_miss5',
                 'scANVI miss 10%':'label_miss/scanvi_miss10',
                 'scANVI miss %':'label_miss/scanvi_miss',
                 'scPoli miss 5%':'label_miss/scpoli_miss5',
                 'scPoli miss 10%':'label_miss/scpoli_miss10',
                 'scPoli miss %':'label_miss/scpoli_miss',
                  
                 'Palette shuffled 5%':'label_shuffled/Palette_shuffled5',
                 'Palette shuffled 10%':'label_shuffled/Palette_shuffled10',
                 'Palette shuffled %':'label_shuffled/Palette_shuffled',
                 'ssSTACAS shuffled 5%':'label_shuffled/ssSTACAS_shuffled5',
                 'ssSTACAS shuffled 10%':'label_shuffled/ssSTACAS_shuffled10',
                 'ssSTACAS shuffled %':'label_shuffled/ssSTACAS_shuffled',
                 'SIGNAL shuffled 5%':'label_shuffled/SIGNAL_shuffled5',
                 'SIGNAL shuffled 10%':'label_shuffled/SIGNAL_shuffled10',
                 'SIGNAL shuffled %':'label_shuffled/SIGNAL_shuffled',
                 'scPoli shuffled 5%':'label_shuffled/scpoli_shuffled5',
                 'scPoli shuffled 10%':'label_shuffled/scpoli_shuffled10',
                 'scPoli shuffled %':'label_shuffled/scpoli_shuffled',
                 'scANVI shuffled 5%':'label_shuffled/scanvi_shuffled5',
                 'scANVI shuffled 10%':'label_shuffled/scanvi_shuffled10',
                 'scANVI shuffled %':'label_shuffled/scanvi_shuffled'
                 }

condition_key = 'Batch'
cell_type_key = 'CellType'

meta = pyreadr.read_r(meta_dir + "meta.rds")[None]
# dataset = dataset[None]
dataset = np.random.rand(meta.shape[0], )
adata = anndata.AnnData(X=dataset, obs=meta)

data_dir = out_dir
for i in method_mapping.keys():
    res = pd.read_csv(data_dir + method_mapping[i] + ".csv")
    print(i, "dimension: ", res.shape[0], res.shape[1])
    adata.obsm[i] = res.values[:,:]

batchcorr = BatchCorrection(pcr_comparison = False)
bioconser = BioConservation(isolated_labels = False)

bm = Benchmarker(
    adata,
    batch_key=condition_key,
    label_key=cell_type_key,
    embedding_obsm_keys=list(method_mapping.keys()),
    bio_conservation_metrics = bioconser,
    batch_correction_metrics = batchcorr,
    n_jobs=6,
    )
bm.benchmark()


# In[7]:


bm.plot_results_table(min_max_scale=False)
df = bm.get_results(min_max_scale=False)
df

latent_data = df
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    latent_df = latent_data

latent_df.to_csv('./metrics.csv', index=True)
