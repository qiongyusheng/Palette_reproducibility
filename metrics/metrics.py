import scanpy as sc
import numpy as np
import pandas as pd
import pyreadr
import anndata
import warnings
warnings.filterwarnings('ignore')
from scib_metrics.benchmark import Benchmarker, BatchCorrection, BioConservation

out_dir = '../'
meta_dir = '../'

method_mapping = { "Palette":"Palette",
                  'Multigrate':'Multigrate',
                  'scMoMaT':'scMoMaT',
                  'midas':'midas',
                  'scVAEIT':'scVAEIT',
                  "StabMap":"stab"}
condition_key = 'Batch'
cell_type_key = 'celltype'

meta = pyreadr.read_r(meta_dir + "meta.rds")
meta = meta[None]
# dataset = dataset[None]
dataset = np.random.rand(meta.shape[0], 50)
adata = anndata.AnnData(X=dataset, obs=meta)

data_dir = out_dir
for i in method_mapping.keys():
    res = pd.read_csv(out_dir + method_mapping[i] + ".csv")
    adata.obsm[i] = res.values[:,:20]

# scIB metrics
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

# modality mixing
condition_key = 'modality'
cell_type_key = 'celltype'

## modality kBET
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels = False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

bm1 = Benchmarker(
    adata,
    batch_key=condition_key,
    label_key=cell_type_key,
    embedding_obsm_keys=list(method_mapping.keys()),
    # pre_integrated_embedding_obsm_key = 'Raw',
    bio_conservation_metrics = bioconser,
    batch_correction_metrics = batchcorr,
    n_jobs=6,
    )
bm1.benchmark()
df1 = bm1.get_results(min_max_scale=False,)

metric_df = df[['KMeans NMI','KMeans ARI','Silhouette label','cLISI','Silhouette batch','iLISI','KBET','Graph connectivity']]
metric_df['KBET(modality)'] = df1['KBET']
metric_df = metric_df.iloc[:-1]
del df1
del bm1

## SAS
import numbers
from typing import Any, Mapping, Optional, TypeVar, Union

import anndata as ad
import h5py
import numpy as np
import scipy.sparse
import sklearn.neighbors

Array = Union[np.ndarray, scipy.sparse.spmatrix]
BackedArray = Union[h5py.Dataset, ad._core.sparse_dataset.SparseDataset]
AnyArray = Union[Array, BackedArray]
ArrayOrScalar = Union[np.ndarray, numbers.Number]
Kws = Optional[Mapping[str, Any]]
RandomState = Optional[Union[np.random.RandomState, int]]

T = TypeVar("T")  # Generic type var

def get_rs(x: RandomState = None) -> np.random.RandomState:
    r"""
    Get random state object

    Parameters
    ----------
    x
        Object that can be converted to a random state object

    Returns
    -------
    rs
        Random state object
    """
    if isinstance(x, int):
        return np.random.RandomState(x)
    if isinstance(x, np.random.RandomState):
        return x
    return np.random

def seurat_alignment_score(
        x: np.ndarray, y: np.ndarray, neighbor_frac: float = 0.01,
        n_repeats: int = 4, random_state: RandomState = None, **kwargs
) -> float:
    r"""
    Seurat alignment score

    Parameters
    ----------
    x
        Coordinates
    y
        Batch labels
    neighbor_frac
        Nearest neighbor fraction
    n_repeats
        Number of subsampling repeats
    random_state
        Random state
    **kwargs
        Additional keyword arguments are passed to
        :class:`sklearn.neighbors.NearestNeighbors`

    Returns
    -------
    sas
        Seurat alignment score
    """
    rs = get_rs(random_state)
    idx_list = [np.where(y == u)[0] for u in np.unique(y)]
    min_size = min(idx.size for idx in idx_list)
    repeat_scores = []
    # ncell = np.sum([idx.size for idx in idx_list])
    # k = max(round(ncell * neighbor_frac), 1)
    k = 30
    # nn = sklearn.neighbors.NearestNeighbors(
            # n_neighbors=k + 1, **kwargs
        # ).fit(x)
    # nni = nn.kneighbors(x, return_distance=False)
    # same_y_hits = (
            # y[nni[:, 1:]] == np.expand_dims(y, axis=1)
        # ).sum(axis=1).mean()
    # repeat_score = (k - same_y_hits) * len(idx_list) / (k * (len(idx_list) - 1))
    
    
    for _ in range(n_repeats):
        subsample_idx = np.concatenate([
            rs.choice(idx, min_size, replace=False)
            for idx in idx_list
        ])
        subsample_x = x[subsample_idx]
        subsample_y = y[subsample_idx]
        # k = max(round(subsample_idx.size * neighbor_frac), 1)
        nn = sklearn.neighbors.NearestNeighbors(
            n_neighbors=k + 1, **kwargs
        ).fit(subsample_x)
        nni = nn.kneighbors(subsample_x, return_distance=False)
        same_y_hits = (
            subsample_y[nni[:, 1:]] == np.expand_dims(subsample_y, axis=1)
        ).sum(axis=1).mean()
        repeat_score = (k - same_y_hits) * len(idx_list) / (k * (len(idx_list) - 1))
        repeat_scores.append(min(repeat_score, 1))  # score may exceed 1, if same_y_hits is lower than expected by chance
    return np.mean(repeat_scores).item()
    # return repeat_score


metric_df['SAS'] = 0
modal = np.array(adata.obs['modality'])
for i in range(metric_df.shape[0]):
    method = metric_df.index[i]
    sas_tmp = seurat_alignment_score(adata.obsm[method],modal)
    metric_df.loc[method,'SAS'] = sas_tmp

r_df = pyreadr.read_r(data_dir + "f1_lisi_score_l1.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['F1 LISI'] = r_df['F1_LISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['F1 LISI'] = r_df2['F1_LISI']

metric_df['Bio conservation'] = (metric_df['KMeans NMI']+metric_df['KMeans ARI']+metric_df['Silhouette label']+metric_df['cLISI'])/4
metric_df['Batch correction'] = (metric_df['Silhouette batch']+metric_df['iLISI']+metric_df['KBET']+metric_df['Graph connectivity'])/4
metric_df['Total'] = metric_df['Bio conservation']*0.6 + metric_df['Batch correction']*0.4
metric_df['Modality mixing'] = (metric_df['KBET(modality)']+metric_df['SAS']+metric_df['F1 LISI'])/3

metric_df.loc['Metric Type'] = ['Bio conservation','Bio conservation','Bio conservation','Bio conservation',
                                'Batch correction','Batch correction','Batch correction','Batch correction',
                                'Modality mixing','Modality mixing','Modality mixing',
                               'Aggregate score','Aggregate score','Aggregate score','Aggregate score']

latent_data = metric_df
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    latent_df = latent_data

latent_df.to_csv('../metrics.csv', index=True)