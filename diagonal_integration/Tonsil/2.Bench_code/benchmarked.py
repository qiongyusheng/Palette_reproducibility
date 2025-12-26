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
method_mapping = {"Palette":"Palette",
                  "StabMap":"stab",
                  "Multigrate":"Multigrate", 
                  "scMoMaT":"scMoMaT", 
                  "midas":"midas",
                  "scVAEIT":'scVAEIT',
                  "bindSC":"bindsc",
                  "uniPort":'uniport',
                  "UINMF":"UINMF",
                  'fastMNN':"fastMNN",
                  "Seurat":"Seurat",
                  "Harmony":"Harmony",
                  'MaxFuse':'maxfuse',
                  'scConfluence':'scCon'}
condition_key = 'modality'
cell_type_key = 'celltype'
meta = pyreadr.read_r(meta_dir + "meta.rds")
meta = meta[None]
# prepare adata object
dataset = np.random.rand(191896, 50)
adata = anndata.AnnData(X=dataset, obs=meta)
data_dir = out_dir
for i in method_mapping.keys():
    res = pd.read_csv(data_dir + method_mapping[i] + ".csv")
    #res = res[None]
    # print(i, "dimension: ", res.shape[0], res.shape[1])
    print(i, "dimension: ", res.shape[0], res.shape[1])
    adata.obsm[i] = res.values
    #adata.obsm[i] = res
batchcorr = BatchCorrection(pcr_comparison = False)
bioconser = BioConservation(isolated_labels=False)

bm = Benchmarker(
    adata,
    batch_key=condition_key,
    label_key=cell_type_key,
    embedding_obsm_keys=list(method_mapping.keys()),
    # pre_integrated_embedding_obsm_key = 'Raw',
    bio_conservation_metrics = bioconser,
    batch_correction_metrics = batchcorr,
    n_jobs=6,
    )
bm.benchmark()

bm.plot_results_table(min_max_scale=False)
df = bm.get_results(min_max_scale=False)
latent_data = df
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    latent_df = latent_data

latent_df.to_csv('./metrics.csv', index=True)
