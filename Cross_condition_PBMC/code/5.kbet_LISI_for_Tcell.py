import scanpy as sc
import numpy as np
import pandas as pd
import pyreadr
import anndata
import warnings
warnings.filterwarnings('ignore')
from scib_metrics.benchmark import Benchmarker, BatchCorrection, BioConservation


out_dir = './T_cell_analysis/'
meta_dir = './T_cell_analysis/'
method_mapping = { "newDGE":"newDGE1", 
                  "DGE":"DGE1", 
                 'origion':'origion',
                 'random':"random1"}
condition_key = 'condition'
cell_type_key = 'celltype'

meta = pyreadr.read_r(meta_dir + "meta.rds")[None]
meta.shape

meta = pyreadr.read_r(meta_dir + "meta.rds")
dataset = np.random.rand(12919, 50)
meta = meta[None]
adata = anndata.AnnData(X=dataset, obs=meta)

data_dir = out_dir
for i in method_mapping.keys():
    res = pd.read_csv(data_dir + method_mapping[i] + ".csv")
    print(i, "dimension: ", res.shape[0], res.shape[1])
    adata.obsm[i] = res.values


batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=True,
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
    # pre_integrated_embedding_obsm_key = 'Raw',
    bio_conservation_metrics = bioconser,
    batch_correction_metrics = batchcorr,
    n_jobs=6,
    )
bm.benchmark()


df = bm.get_results(min_max_scale=False)

metric_df = df[['KBET','iLISI']]
metric_df = metric_df.iloc[:-1]

latent_data = metric_df
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    latent_df = latent_data

latent_df.to_csv('/home/server/sqy/data/mosaic/result/new_condition/kbetLISI_T_cell1.csv', index=True)

