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
method_mapping = { "Palette_homo":"Palette_homo",
                  "Palette_non":"Palette_non",
                 'Multigrate_homo':'Multigrate_homo',
                  'Multigrate_non':'Multigrate_non',
                  'scMoMaT_homo':'scMoMaT_homo',
                  'scMoMaT_non':'scMoMaT_non',
                 'MIDAS_homo':'midas',
                  'MIDAS_non':'midas_non',
                 'scVAEIT_homo':'scVAEIT_homo',
                  'scVAEIT_non':'scVAEIT_non',
                 "StabMap_non":"stab_non",
                 "StabMap_homo":"stab_homo",
                 'scPoli':'scpoli',
                 'SIGNAL':'SIGNAL',
                 'Seurat':'Seurat',
                 'scANVI':'scanvi',
                 'Harmony':'Harmony',
                  'fastMNN':'fastMNN'}

condition_key = 'Batch'
cell_type_key = 'CellType'
meta = pyreadr.read_r(meta_dir + "meta.rds")
meta = meta[None]
# dataset = dataset[None]
dataset = np.random.rand(meta.shape[0], 50)
adata = anndata.AnnData(X=dataset, obs=meta)

data_dir = out_dir
for i in method_mapping.keys():
    res = pd.read_csv(out_dir + method_mapping[i] + ".csv")
    adata.obsm[i] = res.values[:,:20]

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

df = bm.get_results(min_max_scale=False)
latent_data = df

if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    latent_df = latent_data
latent_df.to_csv(out_dir + 'metrics.csv', index=True)
