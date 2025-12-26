######################### Modality mix ################################
import scanpy as sc
import numpy as np
import pandas as pd
import pyreadr
import anndata
import warnings
warnings.filterwarnings('ignore')
from scib_metrics.benchmark import Benchmarker, BatchCorrection, BioConservation

data_dir = './Cross_species/MOp/output/'
meta_dir = './Cross_species/MOp/output/'
out_dir = './Cross_species/MOp/output/'
method_mapping = {"Palette":"Palette",
                  'Multigrate':'Multigrate',
                  'scMoMaT':'scMoMaT',
                  'MIDAS':'midas',
                  'scVAEIT':'scVAEIT',
                  "StabMap":"stab"}

meta = pyreadr.read_r(meta_dir + "meta.rds")
meta = meta[None]
dataset = np.random.rand(meta.shape[0], 50)
adata = anndata.AnnData(X=dataset, obs=meta)

data_dir = out_dir
for i in method_mapping.keys():
    res = pd.read_csv(data_dir + method_mapping[i] + ".csv")
    adata.obsm[i] = res.values[:,:20]

condition_key = 'Modal'
cell_type_key = 'anno'

batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels = False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

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

metric_df = df[['iLISI','KBET']]
metric_df = metric_df.iloc[:-1]


latent_data = metric_df
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    latent_df = latent_data

latent_df.to_csv(out_dir + 'modality_mix.csv', index=True)

######################### Species mix ################################
import scanpy as sc
import numpy as np
import pandas as pd
import pyreadr
import anndata
import warnings
warnings.filterwarnings('ignore')
from scib_metrics.benchmark import Benchmarker, BatchCorrection, BioConservation

data_dir = './Cross_species/MOp/output/'
meta_dir = './Cross_species/MOp/output/'
out_dir = './Cross_species/MOp/output/'
method_mapping = {"Palette":"Palette",
                  'Multigrate':'Multigrate',
                  'scMoMaT':'scMoMaT',
                  'MIDAS':'midas',
                  'scVAEIT':'scVAEIT',
                  "StabMap":"stab"}

meta = pyreadr.read_r(meta_dir + "meta.rds")
meta = meta[None]
dataset = np.random.rand(meta.shape[0], 50)
adata = anndata.AnnData(X=dataset, obs=meta)

data_dir = out_dir
for i in method_mapping.keys():
    res = pd.read_csv(data_dir + method_mapping[i] + ".csv")
    adata.obsm[i] = res.values[:,:20]

condition_key = 'species'
cell_type_key = 'anno'

batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels = False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

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

metric_df = df[['iLISI','KBET']]
metric_df = metric_df.iloc[:-1]


latent_data = metric_df
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    latent_df = latent_data
    
latent_df.to_csv(out_dir + 'species_mix.csv', index=True)


######################### Batch mix ################################
import scanpy as sc
import numpy as np
import pandas as pd
import pyreadr
import anndata
import warnings
warnings.filterwarnings('ignore')
from scib_metrics.benchmark import Benchmarker, BatchCorrection, BioConservation

data_dir = './Cross_species/MOp/output/'
meta_dir = './Cross_species/MOp/output/'
out_dir = './Cross_species/MOp/output/'
method_mapping = {"Palette":"Palette",
                  'Multigrate':'Multigrate',
                  'scMoMaT':'scMoMaT',
                  'MIDAS':'midas',
                  'scVAEIT':'scVAEIT',
                  "StabMap":"stab"}

meta = pyreadr.read_r(meta_dir + "meta.rds")
meta = meta[None]
dataset = np.random.rand(meta.shape[0], 50)
adata = anndata.AnnData(X=dataset, obs=meta)

data_dir = out_dir
for i in method_mapping.keys():
    res = pd.read_csv(data_dir + method_mapping[i] + ".csv")
    adata.obsm[i] = res.values[:,:20]

condition_key = 'Batch'
cell_type_key = 'anno'

batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels = False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

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

metric_df = df[['iLISI','KBET']]
metric_df = metric_df.iloc[:-1]


latent_data = metric_df
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    latent_df = latent_data
    

latent_df.to_csv(out_dir + 'batch_mix.csv', index=True)