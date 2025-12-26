######################### Random 1 ################################
import scanpy as sc
import numpy as np
import pandas as pd
import pyreadr
import anndata
import warnings
warnings.filterwarnings('ignore')
from scib_metrics.benchmark import Benchmarker, BatchCorrection, BioConservation
import sys
sys.path.append('./metrics/')
from metrics_py import seurat_alignment_score

path1 = './Benchmarking/Unsupervised/TEA_s2/output/random1/'
path2 = './Benchmarking/Unsupervised/TEA_s2/output/random1/'

out_dir = './Benchmarking/Unsupervised/TEA_s2/output/random1/'
meta_dir = './Benchmarking/Unsupervised/TEA_s2/output/random1/'

method_mapping = { "Palette":"Palette_v1",
                  'Multigrate':'Multigrate_v1',
                  'scMoMaT':'scMoMaT_v1',
                  'MIDAS':'midas',
                  'scVAEIT':'scVAEIT_v1',
                  "StabMap":"stab_v1"}
condition_key = 'batch'
cell_type_key = 'celltype'

meta = pyreadr.read_r(meta_dir + "meta_v1.rds")
meta = meta[None]
# dataset = dataset[None]
dataset = np.random.rand(meta.shape[0], 50)
adata = anndata.AnnData(X=dataset, obs=meta)

data_dir = out_dir
for i in method_mapping.keys():
    if i in ['Palette']:
        res = pd.read_csv(path1 + method_mapping[i] + ".csv")
    else:
        res = pd.read_csv(path2 + method_mapping[i] + ".csv")
    adata.obsm[i] = res.values

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
metric_df = df[['KMeans NMI','KMeans ARI','Silhouette label','cLISI','Silhouette batch','iLISI','KBET','Graph connectivity']]
metric_df = metric_df.iloc[:-1]

if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
    ct_list = list(adata.obs[cell_type_key].cat.categories)
else:
    ct_list = sorted(adata.obs[cell_type_key].unique().tolist())
# 4) 初始化结果表（method x celltype）
df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)

# 5) 准备 Benchmarker 的 metric 对象
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels=False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

# 6) 对每个细胞类型进行判断并可能运行 benchmark
for ct in ct_list:
    adata_ct = adata[adata.obs[cell_type_key] == ct].copy()
    unique_batches = adata_ct.obs[condition_key].unique()
    n_batches = len(unique_batches)

    if n_batches < 2:
        print(f"celltype {ct}: only {n_batches} batch(es) ({unique_batches}), skip KBET")
        continue

    print(f"celltype {ct}: {n_batches} batches -> running benchmark")
    bm = Benchmarker(
        adata_ct,
        batch_key=condition_key,
        label_key=cell_type_key,
        embedding_obsm_keys=list(method_mapping.keys()),
        bio_conservation_metrics=bioconser,
        batch_correction_metrics=batchcorr,
        n_jobs=6,
    )
    bm.benchmark()
    tmp = bm.get_results(min_max_scale=False)

    for method in method_mapping.keys():
        if method in tmp.index and 'KBET' in tmp.columns:
            df.loc[method, ct] = tmp.loc[method, 'KBET']
        else:
            df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)

metric_df['CkBET'] = row_medians

df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)

for method in df.index:
    if method not in adata.obsm:
        print(f"Warning: method '{method}' not found in adata.obsm, skip.")
        continue

    for ct in ct_list:
        adata_ct = adata[adata.obs[cell_type_key] == ct].copy()
        if adata_ct.n_obs == 0:
            df.loc[method, ct] = np.nan
            continue

        unique_batches = adata_ct.obs[condition_key].unique()
        if len(unique_batches) < 2:
            df.loc[method, ct] = np.nan
            continue

        emb = adata_ct.obsm[method]
        if emb.shape[0] != adata_ct.n_obs:
            if emb.shape[1] == adata_ct.n_obs:
                emb = emb.T
            else:
                raise ValueError(
                    f"Embedding shape {emb.shape} 不匹配 adata_ct.n_obs={adata_ct.n_obs} "
                    f"（method={method}, celltype={ct}）。"
                )

        modal = np.array(adata_ct.obs[condition_key])
        n_samples = adata_ct.n_obs

        if n_samples <= 31:
            df.loc[method, ct] = np.nan
            continue

        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[method, ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                k_try = max(1, min(30, n_samples - 1))
                try:
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[method, ct] = sas_tmp
                    print(f"{method}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{method}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[method, ct] = np.nan
            else:
                print(f"Error computing SAS for method {method}, celltype {ct}: {e}")
                df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CSAS'] = row_medians

r_df = pyreadr.read_r(out_dir + "lisi_score_v1.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['CiLISI'] = r_df['CiLISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['CiLISI'] = r_df2['CiLISI']

metric_df['Bio conservation'] = (metric_df['KMeans NMI']+metric_df['KMeans ARI']+metric_df['Silhouette label']+metric_df['cLISI'])/4
metric_df['Batch correction'] = (metric_df['Silhouette batch']+metric_df['iLISI']+metric_df['KBET']+metric_df['Graph connectivity'])/4
metric_df['Total'] = metric_df['Bio conservation']*0.6 + metric_df['Batch correction']*0.4 
metric_df['Modality mixing'] = (metric_df['CkBET']+metric_df['CSAS']+metric_df['CiLISI'])/3

metric_df.loc['Metric Type'] = ['Bio conservation','Bio conservation','Bio conservation','Bio conservation',
                                'Batch correction','Batch correction','Batch correction','Batch correction',
                                'Modality mixing','Modality mixing','Modality mixing',
                               'Aggregate score','Aggregate score','Aggregate score',
                               'Aggregate score']

latent_data = metric_df

if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    latent_df = latent_data

latent_df.to_csv(out_dir + 'metrics_v1.csv', index=True)

######################### Random 2 ################################
import scanpy as sc
import numpy as np
import pandas as pd
import pyreadr
import anndata
import warnings
warnings.filterwarnings('ignore')
from scib_metrics.benchmark import Benchmarker, BatchCorrection, BioConservation
import sys
sys.path.append('./metrics/')
from metrics_py import seurat_alignment_score

path1 = './Benchmarking/Unsupervised/TEA_s2/output/random3/'
path2 = './Benchmarking/Unsupervised/TEA_s2/output/random3/'

out_dir = './Benchmarking/Unsupervised/TEA_s2/output/random3/'
meta_dir = './Benchmarking/Unsupervised/TEA_s2/output/random3/'

method_mapping = { "Palette":"Palette_v2",
                  'Multigrate':'Multigrate_v2',
                  'scMoMaT':'scMoMaT_v2',
                  'MIDAS':'midas_v2',
                  'scVAEIT':'scVAEIT_v2',
                  "StabMap":"stab_v2"}
condition_key = 'batch'
cell_type_key = 'celltype'

meta = pyreadr.read_r(meta_dir + "meta_v2.rds")
meta = meta[None]
# dataset = dataset[None]
dataset = np.random.rand(meta.shape[0], 50)
adata = anndata.AnnData(X=dataset, obs=meta)

data_dir = out_dir
for i in method_mapping.keys():
    if i in ['Palette']:
        res = pd.read_csv(path1 + method_mapping[i] + ".csv")
    else:
        res = pd.read_csv(path2 + method_mapping[i] + ".csv")
    adata.obsm[i] = res.values

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
metric_df = df[['KMeans NMI','KMeans ARI','Silhouette label','cLISI','Silhouette batch','iLISI','KBET','Graph connectivity']]
metric_df = metric_df.iloc[:-1]

if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
    ct_list = list(adata.obs[cell_type_key].cat.categories)
else:
    ct_list = sorted(adata.obs[cell_type_key].unique().tolist())
# 4) 初始化结果表（method x celltype）
df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)

# 5) 准备 Benchmarker 的 metric 对象
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels=False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

# 6) 对每个细胞类型进行判断并可能运行 benchmark
for ct in ct_list:
    adata_ct = adata[adata.obs[cell_type_key] == ct].copy()
    unique_batches = adata_ct.obs[condition_key].unique()
    n_batches = len(unique_batches)

    if n_batches < 2:
        print(f"celltype {ct}: only {n_batches} batch(es) ({unique_batches}), skip KBET")
        continue

    print(f"celltype {ct}: {n_batches} batches -> running benchmark")
    bm = Benchmarker(
        adata_ct,
        batch_key=condition_key,
        label_key=cell_type_key,
        embedding_obsm_keys=list(method_mapping.keys()),
        bio_conservation_metrics=bioconser,
        batch_correction_metrics=batchcorr,
        n_jobs=6,
    )
    bm.benchmark()
    tmp = bm.get_results(min_max_scale=False)

    for method in method_mapping.keys():
        if method in tmp.index and 'KBET' in tmp.columns:
            df.loc[method, ct] = tmp.loc[method, 'KBET']
        else:
            df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)

metric_df['CkBET'] = row_medians

df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)

for method in df.index:
    if method not in adata.obsm:
        print(f"Warning: method '{method}' not found in adata.obsm, skip.")
        continue

    for ct in ct_list:
        adata_ct = adata[adata.obs[cell_type_key] == ct].copy()
        if adata_ct.n_obs == 0:
            df.loc[method, ct] = np.nan
            continue

        unique_batches = adata_ct.obs[condition_key].unique()
        if len(unique_batches) < 2:
            df.loc[method, ct] = np.nan
            continue

        emb = adata_ct.obsm[method]
        if emb.shape[0] != adata_ct.n_obs:
            if emb.shape[1] == adata_ct.n_obs:
                emb = emb.T
            else:
                raise ValueError(
                    f"Embedding shape {emb.shape} 不匹配 adata_ct.n_obs={adata_ct.n_obs} "
                    f"（method={method}, celltype={ct}）。"
                )

        modal = np.array(adata_ct.obs[condition_key])
        n_samples = adata_ct.n_obs

        if n_samples <= 31:
            df.loc[method, ct] = np.nan
            continue

        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[method, ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                k_try = max(1, min(30, n_samples - 1))
                try:
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[method, ct] = sas_tmp
                    print(f"{method}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{method}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[method, ct] = np.nan
            else:
                print(f"Error computing SAS for method {method}, celltype {ct}: {e}")
                df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CSAS'] = row_medians

r_df = pyreadr.read_r(out_dir + "lisi_score_v2.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['CiLISI'] = r_df['CiLISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['CiLISI'] = r_df2['CiLISI']

metric_df['Bio conservation'] = (metric_df['KMeans NMI']+metric_df['KMeans ARI']+metric_df['Silhouette label']+metric_df['cLISI'])/4
metric_df['Batch correction'] = (metric_df['Silhouette batch']+metric_df['iLISI']+metric_df['KBET']+metric_df['Graph connectivity'])/4
metric_df['Total'] = metric_df['Bio conservation']*0.6 + metric_df['Batch correction']*0.4 
metric_df['Modality mixing'] = (metric_df['CkBET']+metric_df['CSAS']+metric_df['CiLISI'])/3

metric_df.loc['Metric Type'] = ['Bio conservation','Bio conservation','Bio conservation','Bio conservation',
                                'Batch correction','Batch correction','Batch correction','Batch correction',
                                'Modality mixing','Modality mixing','Modality mixing',
                               'Aggregate score','Aggregate score','Aggregate score',
                               'Aggregate score']

latent_data = metric_df

if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    latent_df = latent_data

latent_df.to_csv(out_dir + 'metrics_v2.csv', index=True)

######################### Random 3 ################################
import scanpy as sc
import numpy as np
import pandas as pd
import pyreadr
import anndata
import warnings
warnings.filterwarnings('ignore')
from scib_metrics.benchmark import Benchmarker, BatchCorrection, BioConservation
import sys
sys.path.append('./metrics/')
from metrics_py import seurat_alignment_score

path1 = './Benchmarking/Unsupervised/TEA_s2/output/random2/'
path2 = './Benchmarking/Unsupervised/TEA_s2/output/random2/'

out_dir = './Benchmarking/Unsupervised/TEA_s2/output/random2/'
meta_dir = './Benchmarking/Unsupervised/TEA_s2/output/random2/'

method_mapping = { "Palette":"Palette",
                  'Multigrate':'Multigrate',
                  'scMoMaT':'scMoMaT',
                  'MIDAS':'midas',
                  'scVAEIT':'scVAEIT',
                  "StabMap":"stab"}
condition_key = 'batch'
cell_type_key = 'celltype'

meta = pyreadr.read_r(meta_dir + "meta.rds")
meta = meta[None]
# dataset = dataset[None]
dataset = np.random.rand(meta.shape[0], 50)
adata = anndata.AnnData(X=dataset, obs=meta)

data_dir = out_dir
for i in method_mapping.keys():
    if i in ['Palette']:
        res = pd.read_csv(path1 + method_mapping[i] + ".csv")
    else:
        res = pd.read_csv(path2 + method_mapping[i] + ".csv")
    adata.obsm[i] = res.values

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
metric_df = df[['KMeans NMI','KMeans ARI','Silhouette label','cLISI','Silhouette batch','iLISI','KBET','Graph connectivity']]
metric_df = metric_df.iloc[:-1]

if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
    ct_list = list(adata.obs[cell_type_key].cat.categories)
else:
    ct_list = sorted(adata.obs[cell_type_key].unique().tolist())
# 4) 初始化结果表（method x celltype）
df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)

# 5) 准备 Benchmarker 的 metric 对象
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels=False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

# 6) 对每个细胞类型进行判断并可能运行 benchmark
for ct in ct_list:
    adata_ct = adata[adata.obs[cell_type_key] == ct].copy()
    unique_batches = adata_ct.obs[condition_key].unique()
    n_batches = len(unique_batches)

    if n_batches < 2:
        print(f"celltype {ct}: only {n_batches} batch(es) ({unique_batches}), skip KBET")
        continue

    print(f"celltype {ct}: {n_batches} batches -> running benchmark")
    bm = Benchmarker(
        adata_ct,
        batch_key=condition_key,
        label_key=cell_type_key,
        embedding_obsm_keys=list(method_mapping.keys()),
        bio_conservation_metrics=bioconser,
        batch_correction_metrics=batchcorr,
        n_jobs=6,
    )
    bm.benchmark()
    tmp = bm.get_results(min_max_scale=False)

    for method in method_mapping.keys():
        if method in tmp.index and 'KBET' in tmp.columns:
            df.loc[method, ct] = tmp.loc[method, 'KBET']
        else:
            df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)

metric_df['CkBET'] = row_medians

df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)

for method in df.index:
    if method not in adata.obsm:
        print(f"Warning: method '{method}' not found in adata.obsm, skip.")
        continue

    for ct in ct_list:
        adata_ct = adata[adata.obs[cell_type_key] == ct].copy()
        if adata_ct.n_obs == 0:
            df.loc[method, ct] = np.nan
            continue

        unique_batches = adata_ct.obs[condition_key].unique()
        if len(unique_batches) < 2:
            df.loc[method, ct] = np.nan
            continue

        emb = adata_ct.obsm[method]
        if emb.shape[0] != adata_ct.n_obs:
            if emb.shape[1] == adata_ct.n_obs:
                emb = emb.T
            else:
                raise ValueError(
                    f"Embedding shape {emb.shape} 不匹配 adata_ct.n_obs={adata_ct.n_obs} "
                    f"（method={method}, celltype={ct}）。"
                )

        modal = np.array(adata_ct.obs[condition_key])
        n_samples = adata_ct.n_obs

        if n_samples <= 31:
            df.loc[method, ct] = np.nan
            continue

        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[method, ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                k_try = max(1, min(30, n_samples - 1))
                try:
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[method, ct] = sas_tmp
                    print(f"{method}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{method}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[method, ct] = np.nan
            else:
                print(f"Error computing SAS for method {method}, celltype {ct}: {e}")
                df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CSAS'] = row_medians

r_df = pyreadr.read_r(out_dir + "lisi_score.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['CiLISI'] = r_df['CiLISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['CiLISI'] = r_df2['CiLISI']

metric_df['Bio conservation'] = (metric_df['KMeans NMI']+metric_df['KMeans ARI']+metric_df['Silhouette label']+metric_df['cLISI'])/4
metric_df['Batch correction'] = (metric_df['Silhouette batch']+metric_df['iLISI']+metric_df['KBET']+metric_df['Graph connectivity'])/4
metric_df['Total'] = metric_df['Bio conservation']*0.6 + metric_df['Batch correction']*0.4 
metric_df['Modality mixing'] = (metric_df['CkBET']+metric_df['CSAS']+metric_df['CiLISI'])/3

metric_df.loc['Metric Type'] = ['Bio conservation','Bio conservation','Bio conservation','Bio conservation',
                                'Batch correction','Batch correction','Batch correction','Batch correction',
                                'Modality mixing','Modality mixing','Modality mixing',
                               'Aggregate score','Aggregate score','Aggregate score',
                               'Aggregate score']

latent_data = metric_df

if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    latent_df = latent_data

latent_df.to_csv(out_dir + 'metrics.csv', index=True)

######################### Random 4 ################################
import scanpy as sc
import numpy as np
import pandas as pd
import pyreadr
import anndata
import warnings
warnings.filterwarnings('ignore')
from scib_metrics.benchmark import Benchmarker, BatchCorrection, BioConservation
import sys
sys.path.append('./metrics/')
from metrics_py import seurat_alignment_score

path1 = './Benchmarking/Unsupervised/TEA_s2/output/random4/'
path2 = './Benchmarking/Unsupervised/TEA_s2/output/random4/'

out_dir = './Benchmarking/Unsupervised/TEA_s2/output/random4/'
meta_dir = './Benchmarking/Unsupervised/TEA_s2/output/random4/'

method_mapping = { "Palette":"Palette_v3",
                  'Multigrate':'Multigrate_v3',
                  'scMoMaT':'scMoMaT_v3',
                  'MIDAS':'midas_v3',
                  'scVAEIT':'scVAEIT_v3',
                  "StabMap":"stab_v3"}
condition_key = 'batch'
cell_type_key = 'celltype'

meta = pyreadr.read_r(meta_dir + "meta_v3.rds")
meta = meta[None]
# dataset = dataset[None]
dataset = np.random.rand(meta.shape[0], 50)
adata = anndata.AnnData(X=dataset, obs=meta)

data_dir = out_dir
for i in method_mapping.keys():
    if i in ['Palette']:
        res = pd.read_csv(path1 + method_mapping[i] + ".csv")
    else:
        res = pd.read_csv(path2 + method_mapping[i] + ".csv")
    adata.obsm[i] = res.values

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
metric_df = df[['KMeans NMI','KMeans ARI','Silhouette label','cLISI','Silhouette batch','iLISI','KBET','Graph connectivity']]
metric_df = metric_df.iloc[:-1]

if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
    ct_list = list(adata.obs[cell_type_key].cat.categories)
else:
    ct_list = sorted(adata.obs[cell_type_key].unique().tolist())
# 4) 初始化结果表（method x celltype）
df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)

# 5) 准备 Benchmarker 的 metric 对象
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels=False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

# 6) 对每个细胞类型进行判断并可能运行 benchmark
for ct in ct_list:
    adata_ct = adata[adata.obs[cell_type_key] == ct].copy()
    unique_batches = adata_ct.obs[condition_key].unique()
    n_batches = len(unique_batches)

    if n_batches < 2:
        print(f"celltype {ct}: only {n_batches} batch(es) ({unique_batches}), skip KBET")
        continue

    print(f"celltype {ct}: {n_batches} batches -> running benchmark")
    bm = Benchmarker(
        adata_ct,
        batch_key=condition_key,
        label_key=cell_type_key,
        embedding_obsm_keys=list(method_mapping.keys()),
        bio_conservation_metrics=bioconser,
        batch_correction_metrics=batchcorr,
        n_jobs=6,
    )
    bm.benchmark()
    tmp = bm.get_results(min_max_scale=False)

    for method in method_mapping.keys():
        if method in tmp.index and 'KBET' in tmp.columns:
            df.loc[method, ct] = tmp.loc[method, 'KBET']
        else:
            df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)

metric_df['CkBET'] = row_medians

df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)

for method in df.index:
    if method not in adata.obsm:
        print(f"Warning: method '{method}' not found in adata.obsm, skip.")
        continue

    for ct in ct_list:
        adata_ct = adata[adata.obs[cell_type_key] == ct].copy()
        if adata_ct.n_obs == 0:
            df.loc[method, ct] = np.nan
            continue

        unique_batches = adata_ct.obs[condition_key].unique()
        if len(unique_batches) < 2:
            df.loc[method, ct] = np.nan
            continue

        emb = adata_ct.obsm[method]
        if emb.shape[0] != adata_ct.n_obs:
            if emb.shape[1] == adata_ct.n_obs:
                emb = emb.T
            else:
                raise ValueError(
                    f"Embedding shape {emb.shape} 不匹配 adata_ct.n_obs={adata_ct.n_obs} "
                    f"（method={method}, celltype={ct}）。"
                )

        modal = np.array(adata_ct.obs[condition_key])
        n_samples = adata_ct.n_obs

        if n_samples <= 31:
            df.loc[method, ct] = np.nan
            continue

        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[method, ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                k_try = max(1, min(30, n_samples - 1))
                try:
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[method, ct] = sas_tmp
                    print(f"{method}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{method}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[method, ct] = np.nan
            else:
                print(f"Error computing SAS for method {method}, celltype {ct}: {e}")
                df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CSAS'] = row_medians

r_df = pyreadr.read_r(out_dir + "lisi_score_v3.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['CiLISI'] = r_df['CiLISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['CiLISI'] = r_df2['CiLISI']

metric_df['Bio conservation'] = (metric_df['KMeans NMI']+metric_df['KMeans ARI']+metric_df['Silhouette label']+metric_df['cLISI'])/4
metric_df['Batch correction'] = (metric_df['Silhouette batch']+metric_df['iLISI']+metric_df['KBET']+metric_df['Graph connectivity'])/4
metric_df['Total'] = metric_df['Bio conservation']*0.6 + metric_df['Batch correction']*0.4 
metric_df['Modality mixing'] = (metric_df['CkBET']+metric_df['CSAS']+metric_df['CiLISI'])/3

metric_df.loc['Metric Type'] = ['Bio conservation','Bio conservation','Bio conservation','Bio conservation',
                                'Batch correction','Batch correction','Batch correction','Batch correction',
                                'Modality mixing','Modality mixing','Modality mixing',
                               'Aggregate score','Aggregate score','Aggregate score',
                               'Aggregate score']

latent_data = metric_df

if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    latent_df = latent_data

latent_df.to_csv(out_dir + 'metrics_v3.csv', index=True)


######################### Random 5 ################################
import scanpy as sc
import numpy as np
import pandas as pd
import pyreadr
import anndata
import warnings
warnings.filterwarnings('ignore')
from scib_metrics.benchmark import Benchmarker, BatchCorrection, BioConservation
import sys
sys.path.append('./metrics/')
from metrics_py import seurat_alignment_score

path1 = './Benchmarking/Unsupervised/TEA_s2/output/random5/'
path2 = './Benchmarking/Unsupervised/TEA_s2/output/random5/'

out_dir = './Benchmarking/Unsupervised/TEA_s2/output/random5/'
meta_dir = './Benchmarking/Unsupervised/TEA_s2/output/random5/'

method_mapping = { "Palette":"Palette_v4",
                  'Multigrate':'Multigrate_v4',
                  'scMoMaT':'scMoMaT_v4',
                  'MIDAS':'midas_v4',
                  'scVAEIT':'scVAEIT_v4',
                  "StabMap":"stab_v4"}
condition_key = 'batch'
cell_type_key = 'celltype'

meta = pyreadr.read_r(meta_dir + "meta_v4.rds")
meta = meta[None]
# dataset = dataset[None]
dataset = np.random.rand(meta.shape[0], 50)
adata = anndata.AnnData(X=dataset, obs=meta)

data_dir = out_dir
for i in method_mapping.keys():
    if i in ['Palette']:
        res = pd.read_csv(path1 + method_mapping[i] + ".csv")
    else:
        res = pd.read_csv(path2 + method_mapping[i] + ".csv")
    adata.obsm[i] = res.values

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
metric_df = df[['KMeans NMI','KMeans ARI','Silhouette label','cLISI','Silhouette batch','iLISI','KBET','Graph connectivity']]
metric_df = metric_df.iloc[:-1]


if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
    ct_list = list(adata.obs[cell_type_key].cat.categories)
else:
    ct_list = sorted(adata.obs[cell_type_key].unique().tolist())
# 4) 初始化结果表（method x celltype）
df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)

# 5) 准备 Benchmarker 的 metric 对象
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels=False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

# 6) 对每个细胞类型进行判断并可能运行 benchmark
for ct in ct_list:
    adata_ct = adata[adata.obs[cell_type_key] == ct].copy()
    unique_batches = adata_ct.obs[condition_key].unique()
    n_batches = len(unique_batches)

    if n_batches < 2:
        print(f"celltype {ct}: only {n_batches} batch(es) ({unique_batches}), skip KBET")
        continue

    print(f"celltype {ct}: {n_batches} batches -> running benchmark")
    bm = Benchmarker(
        adata_ct,
        batch_key=condition_key,
        label_key=cell_type_key,
        embedding_obsm_keys=list(method_mapping.keys()),
        bio_conservation_metrics=bioconser,
        batch_correction_metrics=batchcorr,
        n_jobs=6,
    )
    bm.benchmark()
    tmp = bm.get_results(min_max_scale=False)

    for method in method_mapping.keys():
        if method in tmp.index and 'KBET' in tmp.columns:
            df.loc[method, ct] = tmp.loc[method, 'KBET']
        else:
            df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)

metric_df['CkBET'] = row_medians

df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)

for method in df.index:
    if method not in adata.obsm:
        print(f"Warning: method '{method}' not found in adata.obsm, skip.")
        continue

    for ct in ct_list:
        adata_ct = adata[adata.obs[cell_type_key] == ct].copy()
        if adata_ct.n_obs == 0:
            df.loc[method, ct] = np.nan
            continue

        unique_batches = adata_ct.obs[condition_key].unique()
        if len(unique_batches) < 2:
            df.loc[method, ct] = np.nan
            continue

        emb = adata_ct.obsm[method]
        if emb.shape[0] != adata_ct.n_obs:
            if emb.shape[1] == adata_ct.n_obs:
                emb = emb.T
            else:
                raise ValueError(
                    f"Embedding shape {emb.shape} 不匹配 adata_ct.n_obs={adata_ct.n_obs} "
                    f"（method={method}, celltype={ct}）。"
                )

        modal = np.array(adata_ct.obs[condition_key])
        n_samples = adata_ct.n_obs

        if n_samples <= 31:
            df.loc[method, ct] = np.nan
            continue

        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[method, ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                k_try = max(1, min(30, n_samples - 1))
                try:
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[method, ct] = sas_tmp
                    print(f"{method}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{method}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[method, ct] = np.nan
            else:
                print(f"Error computing SAS for method {method}, celltype {ct}: {e}")
                df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CSAS'] = row_medians

r_df = pyreadr.read_r(out_dir + "lisi_score_v4.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['CiLISI'] = r_df['CiLISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['CiLISI'] = r_df2['CiLISI']

metric_df['Bio conservation'] = (metric_df['KMeans NMI']+metric_df['KMeans ARI']+metric_df['Silhouette label']+metric_df['cLISI'])/4
metric_df['Batch correction'] = (metric_df['Silhouette batch']+metric_df['iLISI']+metric_df['KBET']+metric_df['Graph connectivity'])/4
metric_df['Total'] = metric_df['Bio conservation']*0.6 + metric_df['Batch correction']*0.4 
metric_df['Modality mixing'] = (metric_df['CkBET']+metric_df['CSAS']+metric_df['CiLISI'])/3

metric_df.loc['Metric Type'] = ['Bio conservation','Bio conservation','Bio conservation','Bio conservation',
                                'Batch correction','Batch correction','Batch correction','Batch correction',
                                'Modality mixing','Modality mixing','Modality mixing',
                               'Aggregate score','Aggregate score','Aggregate score',
                               'Aggregate score']

latent_data = metric_df

if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    latent_df = latent_data

latent_df.to_csv(out_dir + 'metrics_v4.csv', index=True)