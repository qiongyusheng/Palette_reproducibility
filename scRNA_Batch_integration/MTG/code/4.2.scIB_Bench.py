import scanpy as sc
import numpy as np
import pandas as pd
import pyreadr
import anndata
import warnings
warnings.filterwarnings('ignore')
from scib_metrics.benchmark import Benchmarker, BatchCorrection, BioConservation

out_dir = './'
meta_dir = './'

method_mapping = {'Palette Unsupervised':"Palette_unsupervised",
                 'Palette Supervised':"Palette_supervised",
                 'STACAS':'STACAS',
                 'ssSTACAS':'ssSTACAS',
                 'SIGNAL':'SIGNAL',
                 'scPoli':'scpoli',
                 'scANVI':'scanvi',
                 'fastMNN':'fastMNN',
                 'Harmony':'Harmony',
                 'Seurat':'Seurat',
                 'RAW':'raw',
                  
                 'Palette miss 5%':'label_miss/Palette_miss5',
                 'Palette miss 10%':'label_miss/Palette_miss10',
                 'Palette miss 20%':'label_miss/Palette_miss20',
                 'SIGNAL miss 5%':'label_miss/SIGNAL_miss5',
                 'SIGNAL miss 10%':'label_miss/SIGNAL_miss10',
                 'SIGNAL miss 20%':'label_miss/SIGNAL_miss20',
                 'ssSTACAS miss 5%':'label_miss/ssSTACAS_miss5',
                 'ssSTACAS miss 10%':'label_miss/ssSTACAS_miss10',
                 'ssSTACAS miss 20%':'label_miss/ssSTACAS_miss20',
                 'scANVI miss 5%':'label_miss/scanvi_miss5',
                 'scANVI miss 10%':'label_miss/scanvi_miss10',
                 'scANVI miss 20%':'label_miss/scanvi_miss20',
                 'scPoli miss 5%':'label_miss/scpoli_miss5',
                 'scPoli miss 10%':'label_miss/scpoli_miss10',
                 'scPoli miss 20%':'label_miss/scpoli_miss20',
                  
                 'Palette shuffled 5%':'label_shuffled/Palette_shuffled5',
                 'Palette shuffled 10%':'label_shuffled/Palette_shuffled10',
                 'Palette shuffled 20%':'label_shuffled/Palette_shuffled20',
                 'ssSTACAS shuffled 5%':'label_shuffled/ssSTACAS_shuffled5',
                 'ssSTACAS shuffled 10%':'label_shuffled/ssSTACAS_shuffled10',
                 'ssSTACAS shuffled 20%':'label_shuffled/ssSTACAS_shuffled20',
                 'SIGNAL shuffled 5%':'label_shuffled/SIGNAL_shuffled5',
                 'SIGNAL shuffled 10%':'label_shuffled/SIGNAL_shuffled10',
                 'SIGNAL shuffled 20%':'label_shuffled/SIGNAL_shuffled20',
                 'scPoli shuffled 5%':'label_shuffled/scpoli_shuffled5',
                 'scPoli shuffled 10%':'label_shuffled/scpoli_shuffled10',
                 'scPoli shuffled 20%':'label_shuffled/scpoli_shuffled20',
                 'scANVI shuffled 5%':'label_shuffled/scanvi_shuffled5',
                 'scANVI shuffled 10%':'label_shuffled/scanvi_shuffled10',
                 'scANVI shuffled 20%':'label_shuffled/scanvi_shuffled20'}

condition_key = 'Batch'
cell_type_key = 'SubClass'

meta = pyreadr.read_r(meta_dir + "meta.rds")[None]
dataset = np.random.rand(meta.shape[0], 20)
adata = anndata.AnnData(X=dataset, obs=meta)

data_dir = out_dir
for i in method_mapping.keys():
    res = pd.read_csv(data_dir + method_mapping[i] + ".csv")
    print(i, "dimension: ", res.shape[0], res.shape[1])
    adata.obsm[i] = res.values[:,:20]

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

df = bm.get_results(min_max_scale=False)

latent_data = df
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    latent_df = latent_data

latent_df.to_csv('./metrics.csv', index=True)
