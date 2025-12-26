######################### ADT ################################
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

out_dir = './Benchmarking/Reference-based_integration/output/Result/adt/'
meta_dir = './Benchmarking/Reference-based_integration/output/Result/adt/'

method_mapping = { "MIDAS":"MIDAS",
                  'Multigrate':'Multigrate',
                  'Palette_harmony':'Palette_harmony',
                  'Palette_harmony_supervised':'Palette_harmony_supervised',
                   'Palette_fastMNN':'Palette_fastMNN',
                   'Palette_fastMNN_supervised':'Palette_fastMNN_supervised',
                   "Palette_Seurat":"Palette_Seurat",
                   "Palette_Seurat_supervised":"Palette_Seurat_supervised",
                   'Palette_denovo':'Palette_denovo',
                   'Palette_denovo_supervised':'Palette_denovo_supervised'}

condition_key = 'batch'
cell_type_key = 'celltype.l2'
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

condition_key = 'modality'
cell_type_key = 'celltype.l2'
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels = False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

metric_df = df[['KMeans NMI','KMeans ARI','Silhouette label','cLISI','Silhouette batch','iLISI','KBET','Graph connectivity']]
metric_df = metric_df.iloc[:-1]

if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
    ct_list = list(adata.obs[cell_type_key].cat.categories)
else:
    ct_list = sorted(adata.obs[cell_type_key].unique().tolist())

# 4) 初始化结果表（method x celltype）
df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)


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

    # 兼容性赋值（若 tmp 中缺少 method 或 KBET 列则填 NaN）
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
            # 只有一个 batch，跳过
            df.loc[method, ct] = np.nan
            continue

        # 取 embedding
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

        # 如果样本数太少（例如 <=1），直接跳过
        if n_samples <= 31:
            df.loc[method, ct] = np.nan
            continue

        # 先直接尝试计算 SAS
        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[method, ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            # 如果报错与 n_neighbors 有关，尝试用较小的邻居数重试（若函数支持该参数）
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                # 选择一个安全的 k：至少 1，最多 30（或 n_samples-1）
                k_try = max(1, min(30, n_samples - 1))
                try:
                    # 仅在该函数接受 n_neighbors 参数时有效
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[method, ct] = sas_tmp
                    print(f"{method}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{method}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[method, ct] = np.nan
            else:
                # 其它错误直接记录 NaN 并打印原因
                print(f"Error computing SAS for method {method}, celltype {ct}: {e}")
                df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CSAS'] = row_medians

r_df = pyreadr.read_r(out_dir + "lisi_score_modality.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['CiLISI'] = r_df['CiLISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['CiLISI'] = r_df2['CiLISI']

## RQ
condition_key = 'RQ'
cell_type_key = 'celltype.l2'
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels = False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
    ct_list = list(adata.obs[cell_type_key].cat.categories)
else:
    ct_list = sorted(adata.obs[cell_type_key].unique().tolist())

# 4) 初始化结果表（method x celltype）
df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)

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

    # 兼容性赋值（若 tmp 中缺少 method 或 KBET 列则填 NaN）
    for method in method_mapping.keys():
        if method in tmp.index and 'KBET' in tmp.columns:
            df.loc[method, ct] = tmp.loc[method, 'KBET']
        else:
            df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CkBET(RQ)'] = row_medians

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
            # 只有一个 batch，跳过
            df.loc[method, ct] = np.nan
            continue

        # 取 embedding
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

        # 如果样本数太少（例如 <=1），直接跳过
        if n_samples <= 31:
            df.loc[method, ct] = np.nan
            continue

        # 先直接尝试计算 SAS
        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[method, ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            # 如果报错与 n_neighbors 有关，尝试用较小的邻居数重试（若函数支持该参数）
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                # 选择一个安全的 k：至少 1，最多 30（或 n_samples-1）
                k_try = max(1, min(30, n_samples - 1))
                try:
                    # 仅在该函数接受 n_neighbors 参数时有效
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[method, ct] = sas_tmp
                    print(f"{method}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{method}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[method, ct] = np.nan
            else:
                # 其它错误直接记录 NaN 并打印原因
                print(f"Error computing SAS for method {method}, celltype {ct}: {e}")
                df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CSAS(RQ)'] = row_medians

r_df = pyreadr.read_r(out_dir + "lisi_score_modality.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['CiLISI(RQ)'] = r_df['CiLISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['CiLISI(RQ)'] = r_df2['CiLISI']

metric_df['Bio conservation'] = (metric_df['KMeans NMI']+metric_df['KMeans ARI']+metric_df['Silhouette label']+metric_df['cLISI'])/4
metric_df['Batch correction'] = (metric_df['Silhouette batch']+metric_df['iLISI']+metric_df['KBET']+metric_df['Graph connectivity'])/4
metric_df['Total'] = metric_df['Bio conservation']*0.6 + metric_df['Batch correction']*0.4
metric_df['Modality mixing'] = (metric_df['CkBET']+metric_df['CSAS']+metric_df['CiLISI'])/3  
metric_df['RQ mixing'] = (metric_df['CkBET(RQ)']+metric_df['CSAS(RQ)']+metric_df['CiLISI(RQ)'])/3  

metric_df.loc['Metric Type'] = ['Bio conservation','Bio conservation','Bio conservation','Bio conservation',
                                'Batch correction','Batch correction','Batch correction','Batch correction',
                                'Modality mixing','Modality mixing','Modality mixing',
                                'RQ mixing','RQ mixing','RQ mixing',
                               'Aggregate score','Aggregate score','Aggregate score',
                               'Aggregate score','Aggregate score']

latent_data = metric_df
# 如果 latent_data 是一个 NumPy 数组，将其转换为 DataFrame
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    # 如果 latent_data 已经是 DataFrame，则直接使用
    latent_df = latent_data

# 保存为 CSV 文件
latent_df.to_csv(out_dir + 'metrics_new.csv', index=True)

######################### ADT + ATAC ################################
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

out_dir = './Benchmarking/Reference-based_integration/output/Result/adt_atac/'
meta_dir = './Benchmarking/Reference-based_integration/output/Result/adt_atac/'

method_mapping = { "MIDAS":"MIDAS",
                  'Multigrate':'Multigrate',
                  'Palette_harmony':'Palette_harmony',
                  'Palette_harmony_supervised':'Palette_harmony_supervised',
                   'Palette_fastMNN':'Palette_fastMNN',
                   'Palette_fastMNN_supervised':'Palette_fastMNN_supervised',
                   "Palette_Seurat":"Palette_Seurat",
                   "Palette_Seurat_supervised":"Palette_Seurat_supervised",
                   'Palette_denovo':'Palette_denovo',
                   'Palette_denovo_supervised':'Palette_denovo_supervised'}

condition_key = 'batch'
cell_type_key = 'celltype.l2'
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

condition_key = 'modality'
cell_type_key = 'celltype.l2'
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels = False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

metric_df = df[['KMeans NMI','KMeans ARI','Silhouette label','cLISI','Silhouette batch','iLISI','KBET','Graph connectivity']]
metric_df = metric_df.iloc[:-1]

if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
    ct_list = list(adata.obs[cell_type_key].cat.categories)
else:
    ct_list = sorted(adata.obs[cell_type_key].unique().tolist())

# 4) 初始化结果表（method x celltype）
df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)


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

    # 兼容性赋值（若 tmp 中缺少 method 或 KBET 列则填 NaN）
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
            # 只有一个 batch，跳过
            df.loc[method, ct] = np.nan
            continue

        # 取 embedding
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

        # 如果样本数太少（例如 <=1），直接跳过
        if n_samples <= 31:
            df.loc[method, ct] = np.nan
            continue

        # 先直接尝试计算 SAS
        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[method, ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            # 如果报错与 n_neighbors 有关，尝试用较小的邻居数重试（若函数支持该参数）
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                # 选择一个安全的 k：至少 1，最多 30（或 n_samples-1）
                k_try = max(1, min(30, n_samples - 1))
                try:
                    # 仅在该函数接受 n_neighbors 参数时有效
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[method, ct] = sas_tmp
                    print(f"{method}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{method}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[method, ct] = np.nan
            else:
                # 其它错误直接记录 NaN 并打印原因
                print(f"Error computing SAS for method {method}, celltype {ct}: {e}")
                df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CSAS'] = row_medians

r_df = pyreadr.read_r(out_dir + "lisi_score_modality.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['CiLISI'] = r_df['CiLISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['CiLISI'] = r_df2['CiLISI']

## RQ
condition_key = 'RQ'
cell_type_key = 'celltype.l2'
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels = False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
    ct_list = list(adata.obs[cell_type_key].cat.categories)
else:
    ct_list = sorted(adata.obs[cell_type_key].unique().tolist())

# 4) 初始化结果表（method x celltype）
df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)

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

    # 兼容性赋值（若 tmp 中缺少 method 或 KBET 列则填 NaN）
    for method in method_mapping.keys():
        if method in tmp.index and 'KBET' in tmp.columns:
            df.loc[method, ct] = tmp.loc[method, 'KBET']
        else:
            df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CkBET(RQ)'] = row_medians

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
            # 只有一个 batch，跳过
            df.loc[method, ct] = np.nan
            continue

        # 取 embedding
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

        # 如果样本数太少（例如 <=1），直接跳过
        if n_samples <= 31:
            df.loc[method, ct] = np.nan
            continue

        # 先直接尝试计算 SAS
        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[method, ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            # 如果报错与 n_neighbors 有关，尝试用较小的邻居数重试（若函数支持该参数）
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                # 选择一个安全的 k：至少 1，最多 30（或 n_samples-1）
                k_try = max(1, min(30, n_samples - 1))
                try:
                    # 仅在该函数接受 n_neighbors 参数时有效
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[method, ct] = sas_tmp
                    print(f"{method}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{method}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[method, ct] = np.nan
            else:
                # 其它错误直接记录 NaN 并打印原因
                print(f"Error computing SAS for method {method}, celltype {ct}: {e}")
                df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CSAS(RQ)'] = row_medians

r_df = pyreadr.read_r(out_dir + "lisi_score_modality.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['CiLISI(RQ)'] = r_df['CiLISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['CiLISI(RQ)'] = r_df2['CiLISI']

metric_df['Bio conservation'] = (metric_df['KMeans NMI']+metric_df['KMeans ARI']+metric_df['Silhouette label']+metric_df['cLISI'])/4
metric_df['Batch correction'] = (metric_df['Silhouette batch']+metric_df['iLISI']+metric_df['KBET']+metric_df['Graph connectivity'])/4
metric_df['Total'] = metric_df['Bio conservation']*0.6 + metric_df['Batch correction']*0.4
metric_df['Modality mixing'] = (metric_df['CkBET']+metric_df['CSAS']+metric_df['CiLISI'])/3  
metric_df['RQ mixing'] = (metric_df['CkBET(RQ)']+metric_df['CSAS(RQ)']+metric_df['CiLISI(RQ)'])/3  

metric_df.loc['Metric Type'] = ['Bio conservation','Bio conservation','Bio conservation','Bio conservation',
                                'Batch correction','Batch correction','Batch correction','Batch correction',
                                'Modality mixing','Modality mixing','Modality mixing',
                                'RQ mixing','RQ mixing','RQ mixing',
                               'Aggregate score','Aggregate score','Aggregate score',
                               'Aggregate score','Aggregate score']

latent_data = metric_df
# 如果 latent_data 是一个 NumPy 数组，将其转换为 DataFrame
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    # 如果 latent_data 已经是 DataFrame，则直接使用
    latent_df = latent_data

# 保存为 CSV 文件
latent_df.to_csv(out_dir + 'metrics_new.csv', index=True)

######################### ATAC ################################
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

out_dir = './Benchmarking/Reference-based_integration/output/Result/atac/'
meta_dir = './Benchmarking/Reference-based_integration/output/Result/atac/'

method_mapping = { "MIDAS":"MIDAS",
                  'Multigrate':'Multigrate',
                  'Palette_harmony':'Palette_harmony',
                  'Palette_harmony_supervised':'Palette_harmony_supervised',
                   'Palette_fastMNN':'Palette_fastMNN',
                   'Palette_fastMNN_supervised':'Palette_fastMNN_supervised',
                   "Palette_Seurat":"Palette_Seurat",
                   "Palette_Seurat_supervised":"Palette_Seurat_supervised",
                   'Palette_denovo':'Palette_denovo',
                   'Palette_denovo_supervised':'Palette_denovo_supervised'}

condition_key = 'batch'
cell_type_key = 'celltype.l2'
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

condition_key = 'modality'
cell_type_key = 'celltype.l2'
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels = False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

metric_df = df[['KMeans NMI','KMeans ARI','Silhouette label','cLISI','Silhouette batch','iLISI','KBET','Graph connectivity']]
metric_df = metric_df.iloc[:-1]

if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
    ct_list = list(adata.obs[cell_type_key].cat.categories)
else:
    ct_list = sorted(adata.obs[cell_type_key].unique().tolist())

# 4) 初始化结果表（method x celltype）
df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)


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

    # 兼容性赋值（若 tmp 中缺少 method 或 KBET 列则填 NaN）
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
            # 只有一个 batch，跳过
            df.loc[method, ct] = np.nan
            continue

        # 取 embedding
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

        # 如果样本数太少（例如 <=1），直接跳过
        if n_samples <= 31:
            df.loc[method, ct] = np.nan
            continue

        # 先直接尝试计算 SAS
        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[method, ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            # 如果报错与 n_neighbors 有关，尝试用较小的邻居数重试（若函数支持该参数）
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                # 选择一个安全的 k：至少 1，最多 30（或 n_samples-1）
                k_try = max(1, min(30, n_samples - 1))
                try:
                    # 仅在该函数接受 n_neighbors 参数时有效
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[method, ct] = sas_tmp
                    print(f"{method}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{method}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[method, ct] = np.nan
            else:
                # 其它错误直接记录 NaN 并打印原因
                print(f"Error computing SAS for method {method}, celltype {ct}: {e}")
                df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CSAS'] = row_medians

r_df = pyreadr.read_r(out_dir + "lisi_score_modality.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['CiLISI'] = r_df['CiLISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['CiLISI'] = r_df2['CiLISI']

## RQ
condition_key = 'RQ'
cell_type_key = 'celltype.l2'
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels = False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
    ct_list = list(adata.obs[cell_type_key].cat.categories)
else:
    ct_list = sorted(adata.obs[cell_type_key].unique().tolist())

# 4) 初始化结果表（method x celltype）
df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)

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

    # 兼容性赋值（若 tmp 中缺少 method 或 KBET 列则填 NaN）
    for method in method_mapping.keys():
        if method in tmp.index and 'KBET' in tmp.columns:
            df.loc[method, ct] = tmp.loc[method, 'KBET']
        else:
            df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CkBET(RQ)'] = row_medians

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
            # 只有一个 batch，跳过
            df.loc[method, ct] = np.nan
            continue

        # 取 embedding
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

        # 如果样本数太少（例如 <=1），直接跳过
        if n_samples <= 31:
            df.loc[method, ct] = np.nan
            continue

        # 先直接尝试计算 SAS
        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[method, ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            # 如果报错与 n_neighbors 有关，尝试用较小的邻居数重试（若函数支持该参数）
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                # 选择一个安全的 k：至少 1，最多 30（或 n_samples-1）
                k_try = max(1, min(30, n_samples - 1))
                try:
                    # 仅在该函数接受 n_neighbors 参数时有效
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[method, ct] = sas_tmp
                    print(f"{method}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{method}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[method, ct] = np.nan
            else:
                # 其它错误直接记录 NaN 并打印原因
                print(f"Error computing SAS for method {method}, celltype {ct}: {e}")
                df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CSAS(RQ)'] = row_medians

r_df = pyreadr.read_r(out_dir + "lisi_score_modality.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['CiLISI(RQ)'] = r_df['CiLISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['CiLISI(RQ)'] = r_df2['CiLISI']

metric_df['Bio conservation'] = (metric_df['KMeans NMI']+metric_df['KMeans ARI']+metric_df['Silhouette label']+metric_df['cLISI'])/4
metric_df['Batch correction'] = (metric_df['Silhouette batch']+metric_df['iLISI']+metric_df['KBET']+metric_df['Graph connectivity'])/4
metric_df['Total'] = metric_df['Bio conservation']*0.6 + metric_df['Batch correction']*0.4
metric_df['Modality mixing'] = (metric_df['CkBET']+metric_df['CSAS']+metric_df['CiLISI'])/3  
metric_df['RQ mixing'] = (metric_df['CkBET(RQ)']+metric_df['CSAS(RQ)']+metric_df['CiLISI(RQ)'])/3  

metric_df.loc['Metric Type'] = ['Bio conservation','Bio conservation','Bio conservation','Bio conservation',
                                'Batch correction','Batch correction','Batch correction','Batch correction',
                                'Modality mixing','Modality mixing','Modality mixing',
                                'RQ mixing','RQ mixing','RQ mixing',
                               'Aggregate score','Aggregate score','Aggregate score',
                               'Aggregate score','Aggregate score']

latent_data = metric_df
# 如果 latent_data 是一个 NumPy 数组，将其转换为 DataFrame
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    # 如果 latent_data 已经是 DataFrame，则直接使用
    latent_df = latent_data

# 保存为 CSV 文件
latent_df.to_csv(out_dir + 'metrics_new.csv', index=True)


######################### CITE ################################
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

out_dir = './Benchmarking/Reference-based_integration/output/Result/cite/'
meta_dir = './Benchmarking/Reference-based_integration/output/Result/cite/'

method_mapping = { "MIDAS":"MIDAS",
                  'Multigrate':'Multigrate',
                  'Palette_harmony':'Palette_harmony',
                  'Palette_harmony_supervised':'Palette_harmony_supervised',
                   'Palette_fastMNN':'Palette_fastMNN',
                   'Palette_fastMNN_supervised':'Palette_fastMNN_supervised',
                   "Palette_Seurat":"Palette_Seurat",
                   "Palette_Seurat_supervised":"Palette_Seurat_supervised",
                   'Palette_denovo':'Palette_denovo',
                   'Palette_denovo_supervised':'Palette_denovo_supervised'}

condition_key = 'batch'
cell_type_key = 'celltype.l2'
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

condition_key = 'modality'
cell_type_key = 'celltype.l2'
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels = False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

metric_df = df[['KMeans NMI','KMeans ARI','Silhouette label','cLISI','Silhouette batch','iLISI','KBET','Graph connectivity']]
metric_df = metric_df.iloc[:-1]

if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
    ct_list = list(adata.obs[cell_type_key].cat.categories)
else:
    ct_list = sorted(adata.obs[cell_type_key].unique().tolist())

# 4) 初始化结果表（method x celltype）
df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)


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

    # 兼容性赋值（若 tmp 中缺少 method 或 KBET 列则填 NaN）
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
            # 只有一个 batch，跳过
            df.loc[method, ct] = np.nan
            continue

        # 取 embedding
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

        # 如果样本数太少（例如 <=1），直接跳过
        if n_samples <= 31:
            df.loc[method, ct] = np.nan
            continue

        # 先直接尝试计算 SAS
        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[method, ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            # 如果报错与 n_neighbors 有关，尝试用较小的邻居数重试（若函数支持该参数）
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                # 选择一个安全的 k：至少 1，最多 30（或 n_samples-1）
                k_try = max(1, min(30, n_samples - 1))
                try:
                    # 仅在该函数接受 n_neighbors 参数时有效
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[method, ct] = sas_tmp
                    print(f"{method}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{method}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[method, ct] = np.nan
            else:
                # 其它错误直接记录 NaN 并打印原因
                print(f"Error computing SAS for method {method}, celltype {ct}: {e}")
                df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CSAS'] = row_medians

r_df = pyreadr.read_r(out_dir + "lisi_score_modality.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['CiLISI'] = r_df['CiLISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['CiLISI'] = r_df2['CiLISI']

## RQ
condition_key = 'RQ'
cell_type_key = 'celltype.l2'
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels = False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
    ct_list = list(adata.obs[cell_type_key].cat.categories)
else:
    ct_list = sorted(adata.obs[cell_type_key].unique().tolist())

# 4) 初始化结果表（method x celltype）
df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)

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

    # 兼容性赋值（若 tmp 中缺少 method 或 KBET 列则填 NaN）
    for method in method_mapping.keys():
        if method in tmp.index and 'KBET' in tmp.columns:
            df.loc[method, ct] = tmp.loc[method, 'KBET']
        else:
            df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CkBET(RQ)'] = row_medians

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
            # 只有一个 batch，跳过
            df.loc[method, ct] = np.nan
            continue

        # 取 embedding
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

        # 如果样本数太少（例如 <=1），直接跳过
        if n_samples <= 31:
            df.loc[method, ct] = np.nan
            continue

        # 先直接尝试计算 SAS
        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[method, ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            # 如果报错与 n_neighbors 有关，尝试用较小的邻居数重试（若函数支持该参数）
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                # 选择一个安全的 k：至少 1，最多 30（或 n_samples-1）
                k_try = max(1, min(30, n_samples - 1))
                try:
                    # 仅在该函数接受 n_neighbors 参数时有效
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[method, ct] = sas_tmp
                    print(f"{method}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{method}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[method, ct] = np.nan
            else:
                # 其它错误直接记录 NaN 并打印原因
                print(f"Error computing SAS for method {method}, celltype {ct}: {e}")
                df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CSAS(RQ)'] = row_medians

r_df = pyreadr.read_r(out_dir + "lisi_score_modality.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['CiLISI(RQ)'] = r_df['CiLISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['CiLISI(RQ)'] = r_df2['CiLISI']

metric_df['Bio conservation'] = (metric_df['KMeans NMI']+metric_df['KMeans ARI']+metric_df['Silhouette label']+metric_df['cLISI'])/4
metric_df['Batch correction'] = (metric_df['Silhouette batch']+metric_df['iLISI']+metric_df['KBET']+metric_df['Graph connectivity'])/4
metric_df['Total'] = metric_df['Bio conservation']*0.6 + metric_df['Batch correction']*0.4
metric_df['Modality mixing'] = (metric_df['CkBET']+metric_df['CSAS']+metric_df['CiLISI'])/3  
metric_df['RQ mixing'] = (metric_df['CkBET(RQ)']+metric_df['CSAS(RQ)']+metric_df['CiLISI(RQ)'])/3  

metric_df.loc['Metric Type'] = ['Bio conservation','Bio conservation','Bio conservation','Bio conservation',
                                'Batch correction','Batch correction','Batch correction','Batch correction',
                                'Modality mixing','Modality mixing','Modality mixing',
                                'RQ mixing','RQ mixing','RQ mixing',
                               'Aggregate score','Aggregate score','Aggregate score',
                               'Aggregate score','Aggregate score']

latent_data = metric_df
# 如果 latent_data 是一个 NumPy 数组，将其转换为 DataFrame
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    # 如果 latent_data 已经是 DataFrame，则直接使用
    latent_df = latent_data

# 保存为 CSV 文件
latent_df.to_csv(out_dir + 'metrics_new.csv', index=True)

######################### CITE + Multiome ################################
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

out_dir = './Benchmarking/Reference-based_integration/output/Result/cite_multi/'
meta_dir = './Benchmarking/Reference-based_integration/output/Result/cite_multi/'

method_mapping = { "MIDAS":"MIDAS",
                  'Multigrate':'Multigrate',
                  'Palette_harmony':'Palette_harmony',
                  'Palette_harmony_supervised':'Palette_harmony_supervised',
                   'Palette_fastMNN':'Palette_fastMNN',
                   'Palette_fastMNN_supervised':'Palette_fastMNN_supervised',
                   "Palette_Seurat":"Palette_Seurat",
                   "Palette_Seurat_supervised":"Palette_Seurat_supervised",
                   'Palette_denovo':'Palette_denovo',
                   'Palette_denovo_supervised':'Palette_denovo_supervised'}

condition_key = 'batch'
cell_type_key = 'celltype.l2'
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

condition_key = 'modality'
cell_type_key = 'celltype.l2'
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels = False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

metric_df = df[['KMeans NMI','KMeans ARI','Silhouette label','cLISI','Silhouette batch','iLISI','KBET','Graph connectivity']]
metric_df = metric_df.iloc[:-1]

if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
    ct_list = list(adata.obs[cell_type_key].cat.categories)
else:
    ct_list = sorted(adata.obs[cell_type_key].unique().tolist())

# 4) 初始化结果表（method x celltype）
df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)


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

    # 兼容性赋值（若 tmp 中缺少 method 或 KBET 列则填 NaN）
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
            # 只有一个 batch，跳过
            df.loc[method, ct] = np.nan
            continue

        # 取 embedding
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

        # 如果样本数太少（例如 <=1），直接跳过
        if n_samples <= 31:
            df.loc[method, ct] = np.nan
            continue

        # 先直接尝试计算 SAS
        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[method, ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            # 如果报错与 n_neighbors 有关，尝试用较小的邻居数重试（若函数支持该参数）
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                # 选择一个安全的 k：至少 1，最多 30（或 n_samples-1）
                k_try = max(1, min(30, n_samples - 1))
                try:
                    # 仅在该函数接受 n_neighbors 参数时有效
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[method, ct] = sas_tmp
                    print(f"{method}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{method}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[method, ct] = np.nan
            else:
                # 其它错误直接记录 NaN 并打印原因
                print(f"Error computing SAS for method {method}, celltype {ct}: {e}")
                df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CSAS'] = row_medians

r_df = pyreadr.read_r(out_dir + "lisi_score_modality.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['CiLISI'] = r_df['CiLISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['CiLISI'] = r_df2['CiLISI']

## RQ
condition_key = 'RQ'
cell_type_key = 'celltype.l2'
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels = False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
    ct_list = list(adata.obs[cell_type_key].cat.categories)
else:
    ct_list = sorted(adata.obs[cell_type_key].unique().tolist())

# 4) 初始化结果表（method x celltype）
df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)

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

    # 兼容性赋值（若 tmp 中缺少 method 或 KBET 列则填 NaN）
    for method in method_mapping.keys():
        if method in tmp.index and 'KBET' in tmp.columns:
            df.loc[method, ct] = tmp.loc[method, 'KBET']
        else:
            df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CkBET(RQ)'] = row_medians

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
            # 只有一个 batch，跳过
            df.loc[method, ct] = np.nan
            continue

        # 取 embedding
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

        # 如果样本数太少（例如 <=1），直接跳过
        if n_samples <= 31:
            df.loc[method, ct] = np.nan
            continue

        # 先直接尝试计算 SAS
        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[method, ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            # 如果报错与 n_neighbors 有关，尝试用较小的邻居数重试（若函数支持该参数）
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                # 选择一个安全的 k：至少 1，最多 30（或 n_samples-1）
                k_try = max(1, min(30, n_samples - 1))
                try:
                    # 仅在该函数接受 n_neighbors 参数时有效
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[method, ct] = sas_tmp
                    print(f"{method}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{method}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[method, ct] = np.nan
            else:
                # 其它错误直接记录 NaN 并打印原因
                print(f"Error computing SAS for method {method}, celltype {ct}: {e}")
                df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CSAS(RQ)'] = row_medians

r_df = pyreadr.read_r(out_dir + "lisi_score_modality.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['CiLISI(RQ)'] = r_df['CiLISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['CiLISI(RQ)'] = r_df2['CiLISI']

metric_df['Bio conservation'] = (metric_df['KMeans NMI']+metric_df['KMeans ARI']+metric_df['Silhouette label']+metric_df['cLISI'])/4
metric_df['Batch correction'] = (metric_df['Silhouette batch']+metric_df['iLISI']+metric_df['KBET']+metric_df['Graph connectivity'])/4
metric_df['Total'] = metric_df['Bio conservation']*0.6 + metric_df['Batch correction']*0.4
metric_df['Modality mixing'] = (metric_df['CkBET']+metric_df['CSAS']+metric_df['CiLISI'])/3  
metric_df['RQ mixing'] = (metric_df['CkBET(RQ)']+metric_df['CSAS(RQ)']+metric_df['CiLISI(RQ)'])/3  

metric_df.loc['Metric Type'] = ['Bio conservation','Bio conservation','Bio conservation','Bio conservation',
                                'Batch correction','Batch correction','Batch correction','Batch correction',
                                'Modality mixing','Modality mixing','Modality mixing',
                                'RQ mixing','RQ mixing','RQ mixing',
                               'Aggregate score','Aggregate score','Aggregate score',
                               'Aggregate score','Aggregate score']

latent_data = metric_df
# 如果 latent_data 是一个 NumPy 数组，将其转换为 DataFrame
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    # 如果 latent_data 已经是 DataFrame，则直接使用
    latent_df = latent_data

# 保存为 CSV 文件
latent_df.to_csv(out_dir + 'metrics_new.csv', index=True)


######################### Multiome ################################
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

out_dir = './Benchmarking/Reference-based_integration/output/Result/multiome/'
meta_dir = './Benchmarking/Reference-based_integration/output/Result/multiome/'

method_mapping = { "MIDAS":"MIDAS",
                  'Multigrate':'Multigrate',
                  'Palette_harmony':'Palette_harmony',
                  'Palette_harmony_supervised':'Palette_harmony_supervised',
                   'Palette_fastMNN':'Palette_fastMNN',
                   'Palette_fastMNN_supervised':'Palette_fastMNN_supervised',
                   "Palette_Seurat":"Palette_Seurat",
                   "Palette_Seurat_supervised":"Palette_Seurat_supervised",
                   'Palette_denovo':'Palette_denovo',
                   'Palette_denovo_supervised':'Palette_denovo_supervised'}

condition_key = 'batch'
cell_type_key = 'celltype.l2'
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

condition_key = 'modality'
cell_type_key = 'celltype.l2'
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels = False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

metric_df = df[['KMeans NMI','KMeans ARI','Silhouette label','cLISI','Silhouette batch','iLISI','KBET','Graph connectivity']]
metric_df = metric_df.iloc[:-1]

if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
    ct_list = list(adata.obs[cell_type_key].cat.categories)
else:
    ct_list = sorted(adata.obs[cell_type_key].unique().tolist())

# 4) 初始化结果表（method x celltype）
df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)


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

    # 兼容性赋值（若 tmp 中缺少 method 或 KBET 列则填 NaN）
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
            # 只有一个 batch，跳过
            df.loc[method, ct] = np.nan
            continue

        # 取 embedding
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

        # 如果样本数太少（例如 <=1），直接跳过
        if n_samples <= 31:
            df.loc[method, ct] = np.nan
            continue

        # 先直接尝试计算 SAS
        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[method, ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            # 如果报错与 n_neighbors 有关，尝试用较小的邻居数重试（若函数支持该参数）
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                # 选择一个安全的 k：至少 1，最多 30（或 n_samples-1）
                k_try = max(1, min(30, n_samples - 1))
                try:
                    # 仅在该函数接受 n_neighbors 参数时有效
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[method, ct] = sas_tmp
                    print(f"{method}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{method}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[method, ct] = np.nan
            else:
                # 其它错误直接记录 NaN 并打印原因
                print(f"Error computing SAS for method {method}, celltype {ct}: {e}")
                df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CSAS'] = row_medians

r_df = pyreadr.read_r(out_dir + "lisi_score_modality.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['CiLISI'] = r_df['CiLISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['CiLISI'] = r_df2['CiLISI']

## RQ
condition_key = 'RQ'
cell_type_key = 'celltype.l2'
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels = False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
    ct_list = list(adata.obs[cell_type_key].cat.categories)
else:
    ct_list = sorted(adata.obs[cell_type_key].unique().tolist())

# 4) 初始化结果表（method x celltype）
df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)

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

    # 兼容性赋值（若 tmp 中缺少 method 或 KBET 列则填 NaN）
    for method in method_mapping.keys():
        if method in tmp.index and 'KBET' in tmp.columns:
            df.loc[method, ct] = tmp.loc[method, 'KBET']
        else:
            df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CkBET(RQ)'] = row_medians

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
            # 只有一个 batch，跳过
            df.loc[method, ct] = np.nan
            continue

        # 取 embedding
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

        # 如果样本数太少（例如 <=1），直接跳过
        if n_samples <= 31:
            df.loc[method, ct] = np.nan
            continue

        # 先直接尝试计算 SAS
        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[method, ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            # 如果报错与 n_neighbors 有关，尝试用较小的邻居数重试（若函数支持该参数）
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                # 选择一个安全的 k：至少 1，最多 30（或 n_samples-1）
                k_try = max(1, min(30, n_samples - 1))
                try:
                    # 仅在该函数接受 n_neighbors 参数时有效
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[method, ct] = sas_tmp
                    print(f"{method}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{method}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[method, ct] = np.nan
            else:
                # 其它错误直接记录 NaN 并打印原因
                print(f"Error computing SAS for method {method}, celltype {ct}: {e}")
                df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CSAS(RQ)'] = row_medians

r_df = pyreadr.read_r(out_dir + "lisi_score_modality.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['CiLISI(RQ)'] = r_df['CiLISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['CiLISI(RQ)'] = r_df2['CiLISI']

metric_df['Bio conservation'] = (metric_df['KMeans NMI']+metric_df['KMeans ARI']+metric_df['Silhouette label']+metric_df['cLISI'])/4
metric_df['Batch correction'] = (metric_df['Silhouette batch']+metric_df['iLISI']+metric_df['KBET']+metric_df['Graph connectivity'])/4
metric_df['Total'] = metric_df['Bio conservation']*0.6 + metric_df['Batch correction']*0.4
metric_df['Modality mixing'] = (metric_df['CkBET']+metric_df['CSAS']+metric_df['CiLISI'])/3  
metric_df['RQ mixing'] = (metric_df['CkBET(RQ)']+metric_df['CSAS(RQ)']+metric_df['CiLISI(RQ)'])/3  

metric_df.loc['Metric Type'] = ['Bio conservation','Bio conservation','Bio conservation','Bio conservation',
                                'Batch correction','Batch correction','Batch correction','Batch correction',
                                'Modality mixing','Modality mixing','Modality mixing',
                                'RQ mixing','RQ mixing','RQ mixing',
                               'Aggregate score','Aggregate score','Aggregate score',
                               'Aggregate score','Aggregate score']

latent_data = metric_df
# 如果 latent_data 是一个 NumPy 数组，将其转换为 DataFrame
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    # 如果 latent_data 已经是 DataFrame，则直接使用
    latent_df = latent_data

# 保存为 CSV 文件
latent_df.to_csv(out_dir + 'metrics_new.csv', index=True)


######################### RNA ################################
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

out_dir = './Benchmarking/Reference-based_integration/output/Result/rna/'
meta_dir = './Benchmarking/Reference-based_integration/output/Result/rna/'

method_mapping = { "MIDAS":"MIDAS",
                  'Multigrate':'Multigrate',
                  'Palette_harmony':'Palette_harmony',
                  'Palette_harmony_supervised':'Palette_harmony_supervised',
                   'Palette_fastMNN':'Palette_fastMNN',
                   'Palette_fastMNN_supervised':'Palette_fastMNN_supervised',
                   "Palette_Seurat":"Palette_Seurat",
                   "Palette_Seurat_supervised":"Palette_Seurat_supervised",
                   'Palette_denovo':'Palette_denovo',
                   'Palette_denovo_supervised':'Palette_denovo_supervised'}

condition_key = 'batch'
cell_type_key = 'celltype.l2'
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

condition_key = 'modality'
cell_type_key = 'celltype.l2'
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels = False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

metric_df = df[['KMeans NMI','KMeans ARI','Silhouette label','cLISI','Silhouette batch','iLISI','KBET','Graph connectivity']]
metric_df = metric_df.iloc[:-1]

if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
    ct_list = list(adata.obs[cell_type_key].cat.categories)
else:
    ct_list = sorted(adata.obs[cell_type_key].unique().tolist())

# 4) 初始化结果表（method x celltype）
df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)


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

    # 兼容性赋值（若 tmp 中缺少 method 或 KBET 列则填 NaN）
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
            # 只有一个 batch，跳过
            df.loc[method, ct] = np.nan
            continue

        # 取 embedding
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

        # 如果样本数太少（例如 <=1），直接跳过
        if n_samples <= 31:
            df.loc[method, ct] = np.nan
            continue

        # 先直接尝试计算 SAS
        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[method, ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            # 如果报错与 n_neighbors 有关，尝试用较小的邻居数重试（若函数支持该参数）
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                # 选择一个安全的 k：至少 1，最多 30（或 n_samples-1）
                k_try = max(1, min(30, n_samples - 1))
                try:
                    # 仅在该函数接受 n_neighbors 参数时有效
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[method, ct] = sas_tmp
                    print(f"{method}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{method}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[method, ct] = np.nan
            else:
                # 其它错误直接记录 NaN 并打印原因
                print(f"Error computing SAS for method {method}, celltype {ct}: {e}")
                df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CSAS'] = row_medians

r_df = pyreadr.read_r(out_dir + "lisi_score_modality.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['CiLISI'] = r_df['CiLISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['CiLISI'] = r_df2['CiLISI']

## RQ
condition_key = 'RQ'
cell_type_key = 'celltype.l2'
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels = False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
    ct_list = list(adata.obs[cell_type_key].cat.categories)
else:
    ct_list = sorted(adata.obs[cell_type_key].unique().tolist())

# 4) 初始化结果表（method x celltype）
df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)

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

    # 兼容性赋值（若 tmp 中缺少 method 或 KBET 列则填 NaN）
    for method in method_mapping.keys():
        if method in tmp.index and 'KBET' in tmp.columns:
            df.loc[method, ct] = tmp.loc[method, 'KBET']
        else:
            df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CkBET(RQ)'] = row_medians

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
            # 只有一个 batch，跳过
            df.loc[method, ct] = np.nan
            continue

        # 取 embedding
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

        # 如果样本数太少（例如 <=1），直接跳过
        if n_samples <= 31:
            df.loc[method, ct] = np.nan
            continue

        # 先直接尝试计算 SAS
        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[method, ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            # 如果报错与 n_neighbors 有关，尝试用较小的邻居数重试（若函数支持该参数）
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                # 选择一个安全的 k：至少 1，最多 30（或 n_samples-1）
                k_try = max(1, min(30, n_samples - 1))
                try:
                    # 仅在该函数接受 n_neighbors 参数时有效
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[method, ct] = sas_tmp
                    print(f"{method}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{method}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[method, ct] = np.nan
            else:
                # 其它错误直接记录 NaN 并打印原因
                print(f"Error computing SAS for method {method}, celltype {ct}: {e}")
                df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CSAS(RQ)'] = row_medians

r_df = pyreadr.read_r(out_dir + "lisi_score_modality.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['CiLISI(RQ)'] = r_df['CiLISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['CiLISI(RQ)'] = r_df2['CiLISI']

metric_df['Bio conservation'] = (metric_df['KMeans NMI']+metric_df['KMeans ARI']+metric_df['Silhouette label']+metric_df['cLISI'])/4
metric_df['Batch correction'] = (metric_df['Silhouette batch']+metric_df['iLISI']+metric_df['KBET']+metric_df['Graph connectivity'])/4
metric_df['Total'] = metric_df['Bio conservation']*0.6 + metric_df['Batch correction']*0.4
metric_df['Modality mixing'] = (metric_df['CkBET']+metric_df['CSAS']+metric_df['CiLISI'])/3  
metric_df['RQ mixing'] = (metric_df['CkBET(RQ)']+metric_df['CSAS(RQ)']+metric_df['CiLISI(RQ)'])/3  

metric_df.loc['Metric Type'] = ['Bio conservation','Bio conservation','Bio conservation','Bio conservation',
                                'Batch correction','Batch correction','Batch correction','Batch correction',
                                'Modality mixing','Modality mixing','Modality mixing',
                                'RQ mixing','RQ mixing','RQ mixing',
                               'Aggregate score','Aggregate score','Aggregate score',
                               'Aggregate score','Aggregate score']

latent_data = metric_df
# 如果 latent_data 是一个 NumPy 数组，将其转换为 DataFrame
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    # 如果 latent_data 已经是 DataFrame，则直接使用
    latent_df = latent_data

# 保存为 CSV 文件
latent_df.to_csv(out_dir + 'metrics_new.csv', index=True)

######################### RNA + ADT ################################
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

out_dir = './Benchmarking/Reference-based_integration/output/Result/rna_adt/'
meta_dir = './Benchmarking/Reference-based_integration/output/Result/rna_adt/'

method_mapping = { "MIDAS":"MIDAS",
                  'Multigrate':'Multigrate',
                  'Palette_harmony':'Palette_harmony',
                  'Palette_harmony_supervised':'Palette_harmony_supervised',
                   'Palette_fastMNN':'Palette_fastMNN',
                   'Palette_fastMNN_supervised':'Palette_fastMNN_supervised',
                   "Palette_Seurat":"Palette_Seurat",
                   "Palette_Seurat_supervised":"Palette_Seurat_supervised",
                   'Palette_denovo':'Palette_denovo',
                   'Palette_denovo_supervised':'Palette_denovo_supervised'}

condition_key = 'batch'
cell_type_key = 'celltype.l2'
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

condition_key = 'modality'
cell_type_key = 'celltype.l2'
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels = False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

metric_df = df[['KMeans NMI','KMeans ARI','Silhouette label','cLISI','Silhouette batch','iLISI','KBET','Graph connectivity']]
metric_df = metric_df.iloc[:-1]

if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
    ct_list = list(adata.obs[cell_type_key].cat.categories)
else:
    ct_list = sorted(adata.obs[cell_type_key].unique().tolist())

# 4) 初始化结果表（method x celltype）
df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)


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

    # 兼容性赋值（若 tmp 中缺少 method 或 KBET 列则填 NaN）
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
            # 只有一个 batch，跳过
            df.loc[method, ct] = np.nan
            continue

        # 取 embedding
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

        # 如果样本数太少（例如 <=1），直接跳过
        if n_samples <= 31:
            df.loc[method, ct] = np.nan
            continue

        # 先直接尝试计算 SAS
        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[method, ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            # 如果报错与 n_neighbors 有关，尝试用较小的邻居数重试（若函数支持该参数）
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                # 选择一个安全的 k：至少 1，最多 30（或 n_samples-1）
                k_try = max(1, min(30, n_samples - 1))
                try:
                    # 仅在该函数接受 n_neighbors 参数时有效
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[method, ct] = sas_tmp
                    print(f"{method}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{method}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[method, ct] = np.nan
            else:
                # 其它错误直接记录 NaN 并打印原因
                print(f"Error computing SAS for method {method}, celltype {ct}: {e}")
                df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CSAS'] = row_medians

r_df = pyreadr.read_r(out_dir + "lisi_score_modality.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['CiLISI'] = r_df['CiLISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['CiLISI'] = r_df2['CiLISI']

## RQ
condition_key = 'RQ'
cell_type_key = 'celltype.l2'
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels = False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
    ct_list = list(adata.obs[cell_type_key].cat.categories)
else:
    ct_list = sorted(adata.obs[cell_type_key].unique().tolist())

# 4) 初始化结果表（method x celltype）
df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)

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

    # 兼容性赋值（若 tmp 中缺少 method 或 KBET 列则填 NaN）
    for method in method_mapping.keys():
        if method in tmp.index and 'KBET' in tmp.columns:
            df.loc[method, ct] = tmp.loc[method, 'KBET']
        else:
            df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CkBET(RQ)'] = row_medians

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
            # 只有一个 batch，跳过
            df.loc[method, ct] = np.nan
            continue

        # 取 embedding
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

        # 如果样本数太少（例如 <=1），直接跳过
        if n_samples <= 31:
            df.loc[method, ct] = np.nan
            continue

        # 先直接尝试计算 SAS
        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[method, ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            # 如果报错与 n_neighbors 有关，尝试用较小的邻居数重试（若函数支持该参数）
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                # 选择一个安全的 k：至少 1，最多 30（或 n_samples-1）
                k_try = max(1, min(30, n_samples - 1))
                try:
                    # 仅在该函数接受 n_neighbors 参数时有效
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[method, ct] = sas_tmp
                    print(f"{method}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{method}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[method, ct] = np.nan
            else:
                # 其它错误直接记录 NaN 并打印原因
                print(f"Error computing SAS for method {method}, celltype {ct}: {e}")
                df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CSAS(RQ)'] = row_medians

r_df = pyreadr.read_r(out_dir + "lisi_score_modality.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['CiLISI(RQ)'] = r_df['CiLISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['CiLISI(RQ)'] = r_df2['CiLISI']

metric_df['Bio conservation'] = (metric_df['KMeans NMI']+metric_df['KMeans ARI']+metric_df['Silhouette label']+metric_df['cLISI'])/4
metric_df['Batch correction'] = (metric_df['Silhouette batch']+metric_df['iLISI']+metric_df['KBET']+metric_df['Graph connectivity'])/4
metric_df['Total'] = metric_df['Bio conservation']*0.6 + metric_df['Batch correction']*0.4
metric_df['Modality mixing'] = (metric_df['CkBET']+metric_df['CSAS']+metric_df['CiLISI'])/3  
metric_df['RQ mixing'] = (metric_df['CkBET(RQ)']+metric_df['CSAS(RQ)']+metric_df['CiLISI(RQ)'])/3  

metric_df.loc['Metric Type'] = ['Bio conservation','Bio conservation','Bio conservation','Bio conservation',
                                'Batch correction','Batch correction','Batch correction','Batch correction',
                                'Modality mixing','Modality mixing','Modality mixing',
                                'RQ mixing','RQ mixing','RQ mixing',
                               'Aggregate score','Aggregate score','Aggregate score',
                               'Aggregate score','Aggregate score']

latent_data = metric_df
# 如果 latent_data 是一个 NumPy 数组，将其转换为 DataFrame
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    # 如果 latent_data 已经是 DataFrame，则直接使用
    latent_df = latent_data

# 保存为 CSV 文件
latent_df.to_csv(out_dir + 'metrics_new.csv', index=True)

######################### RNA + ADT + ATAC ################################
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

out_dir = './Benchmarking/Reference-based_integration/output/Result/rna_adt_atac/'
meta_dir = './Benchmarking/Reference-based_integration/output/Result/rna_adt_atac/'

method_mapping = { "MIDAS":"MIDAS",
                  'Multigrate':'Multigrate',
                  'Palette_harmony':'Palette_harmony',
                  'Palette_harmony_supervised':'Palette_harmony_supervised',
                   'Palette_fastMNN':'Palette_fastMNN',
                   'Palette_fastMNN_supervised':'Palette_fastMNN_supervised',
                   "Palette_Seurat":"Palette_Seurat",
                   "Palette_Seurat_supervised":"Palette_Seurat_supervised",
                   'Palette_denovo':'Palette_denovo',
                   'Palette_denovo_supervised':'Palette_denovo_supervised'}

condition_key = 'batch'
cell_type_key = 'celltype.l2'
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

condition_key = 'modality'
cell_type_key = 'celltype.l2'
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels = False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

metric_df = df[['KMeans NMI','KMeans ARI','Silhouette label','cLISI','Silhouette batch','iLISI','KBET','Graph connectivity']]
metric_df = metric_df.iloc[:-1]

if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
    ct_list = list(adata.obs[cell_type_key].cat.categories)
else:
    ct_list = sorted(adata.obs[cell_type_key].unique().tolist())

# 4) 初始化结果表（method x celltype）
df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)


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

    # 兼容性赋值（若 tmp 中缺少 method 或 KBET 列则填 NaN）
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
            # 只有一个 batch，跳过
            df.loc[method, ct] = np.nan
            continue

        # 取 embedding
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

        # 如果样本数太少（例如 <=1），直接跳过
        if n_samples <= 31:
            df.loc[method, ct] = np.nan
            continue

        # 先直接尝试计算 SAS
        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[method, ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            # 如果报错与 n_neighbors 有关，尝试用较小的邻居数重试（若函数支持该参数）
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                # 选择一个安全的 k：至少 1，最多 30（或 n_samples-1）
                k_try = max(1, min(30, n_samples - 1))
                try:
                    # 仅在该函数接受 n_neighbors 参数时有效
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[method, ct] = sas_tmp
                    print(f"{method}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{method}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[method, ct] = np.nan
            else:
                # 其它错误直接记录 NaN 并打印原因
                print(f"Error computing SAS for method {method}, celltype {ct}: {e}")
                df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CSAS'] = row_medians

r_df = pyreadr.read_r(out_dir + "lisi_score_modality.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['CiLISI'] = r_df['CiLISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['CiLISI'] = r_df2['CiLISI']

## RQ
condition_key = 'RQ'
cell_type_key = 'celltype.l2'
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels = False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
    ct_list = list(adata.obs[cell_type_key].cat.categories)
else:
    ct_list = sorted(adata.obs[cell_type_key].unique().tolist())

# 4) 初始化结果表（method x celltype）
df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)

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

    # 兼容性赋值（若 tmp 中缺少 method 或 KBET 列则填 NaN）
    for method in method_mapping.keys():
        if method in tmp.index and 'KBET' in tmp.columns:
            df.loc[method, ct] = tmp.loc[method, 'KBET']
        else:
            df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CkBET(RQ)'] = row_medians

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
            # 只有一个 batch，跳过
            df.loc[method, ct] = np.nan
            continue

        # 取 embedding
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

        # 如果样本数太少（例如 <=1），直接跳过
        if n_samples <= 31:
            df.loc[method, ct] = np.nan
            continue

        # 先直接尝试计算 SAS
        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[method, ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            # 如果报错与 n_neighbors 有关，尝试用较小的邻居数重试（若函数支持该参数）
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                # 选择一个安全的 k：至少 1，最多 30（或 n_samples-1）
                k_try = max(1, min(30, n_samples - 1))
                try:
                    # 仅在该函数接受 n_neighbors 参数时有效
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[method, ct] = sas_tmp
                    print(f"{method}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{method}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[method, ct] = np.nan
            else:
                # 其它错误直接记录 NaN 并打印原因
                print(f"Error computing SAS for method {method}, celltype {ct}: {e}")
                df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CSAS(RQ)'] = row_medians

r_df = pyreadr.read_r(out_dir + "lisi_score_modality.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['CiLISI(RQ)'] = r_df['CiLISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['CiLISI(RQ)'] = r_df2['CiLISI']

metric_df['Bio conservation'] = (metric_df['KMeans NMI']+metric_df['KMeans ARI']+metric_df['Silhouette label']+metric_df['cLISI'])/4
metric_df['Batch correction'] = (metric_df['Silhouette batch']+metric_df['iLISI']+metric_df['KBET']+metric_df['Graph connectivity'])/4
metric_df['Total'] = metric_df['Bio conservation']*0.6 + metric_df['Batch correction']*0.4
metric_df['Modality mixing'] = (metric_df['CkBET']+metric_df['CSAS']+metric_df['CiLISI'])/3  
metric_df['RQ mixing'] = (metric_df['CkBET(RQ)']+metric_df['CSAS(RQ)']+metric_df['CiLISI(RQ)'])/3  

metric_df.loc['Metric Type'] = ['Bio conservation','Bio conservation','Bio conservation','Bio conservation',
                                'Batch correction','Batch correction','Batch correction','Batch correction',
                                'Modality mixing','Modality mixing','Modality mixing',
                                'RQ mixing','RQ mixing','RQ mixing',
                               'Aggregate score','Aggregate score','Aggregate score',
                               'Aggregate score','Aggregate score']

latent_data = metric_df
# 如果 latent_data 是一个 NumPy 数组，将其转换为 DataFrame
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    # 如果 latent_data 已经是 DataFrame，则直接使用
    latent_df = latent_data

# 保存为 CSV 文件
latent_df.to_csv(out_dir + 'metrics_new.csv', index=True)

######################### RNA + ATAC ################################
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

out_dir = './Benchmarking/Reference-based_integration/output/Result/rna_atac/'
meta_dir = './Benchmarking/Reference-based_integration/output/Result/rna_atac/'

method_mapping = { "MIDAS":"MIDAS",
                  'Multigrate':'Multigrate',
                  'Palette_harmony':'Palette_harmony',
                  'Palette_harmony_supervised':'Palette_harmony_supervised',
                   'Palette_fastMNN':'Palette_fastMNN',
                   'Palette_fastMNN_supervised':'Palette_fastMNN_supervised',
                   "Palette_Seurat":"Palette_Seurat",
                   "Palette_Seurat_supervised":"Palette_Seurat_supervised",
                   'Palette_denovo':'Palette_denovo',
                   'Palette_denovo_supervised':'Palette_denovo_supervised'}

condition_key = 'batch'
cell_type_key = 'celltype.l2'
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

condition_key = 'modality'
cell_type_key = 'celltype.l2'
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels = False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

metric_df = df[['KMeans NMI','KMeans ARI','Silhouette label','cLISI','Silhouette batch','iLISI','KBET','Graph connectivity']]
metric_df = metric_df.iloc[:-1]

if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
    ct_list = list(adata.obs[cell_type_key].cat.categories)
else:
    ct_list = sorted(adata.obs[cell_type_key].unique().tolist())

# 4) 初始化结果表（method x celltype）
df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)


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

    # 兼容性赋值（若 tmp 中缺少 method 或 KBET 列则填 NaN）
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
            # 只有一个 batch，跳过
            df.loc[method, ct] = np.nan
            continue

        # 取 embedding
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

        # 如果样本数太少（例如 <=1），直接跳过
        if n_samples <= 31:
            df.loc[method, ct] = np.nan
            continue

        # 先直接尝试计算 SAS
        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[method, ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            # 如果报错与 n_neighbors 有关，尝试用较小的邻居数重试（若函数支持该参数）
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                # 选择一个安全的 k：至少 1，最多 30（或 n_samples-1）
                k_try = max(1, min(30, n_samples - 1))
                try:
                    # 仅在该函数接受 n_neighbors 参数时有效
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[method, ct] = sas_tmp
                    print(f"{method}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{method}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[method, ct] = np.nan
            else:
                # 其它错误直接记录 NaN 并打印原因
                print(f"Error computing SAS for method {method}, celltype {ct}: {e}")
                df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CSAS'] = row_medians

r_df = pyreadr.read_r(out_dir + "lisi_score_modality.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['CiLISI'] = r_df['CiLISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['CiLISI'] = r_df2['CiLISI']

## RQ
condition_key = 'RQ'
cell_type_key = 'celltype.l2'
batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
bioconser = BioConservation(isolated_labels = False,
                            nmi_ari_cluster_labels_kmeans=False,
                            nmi_ari_cluster_labels_leiden=False,
                            silhouette_label=False)

if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
    ct_list = list(adata.obs[cell_type_key].cat.categories)
else:
    ct_list = sorted(adata.obs[cell_type_key].unique().tolist())

# 4) 初始化结果表（method x celltype）
df = pd.DataFrame(np.nan, index=list(method_mapping.keys()), columns=ct_list)

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

    # 兼容性赋值（若 tmp 中缺少 method 或 KBET 列则填 NaN）
    for method in method_mapping.keys():
        if method in tmp.index and 'KBET' in tmp.columns:
            df.loc[method, ct] = tmp.loc[method, 'KBET']
        else:
            df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CkBET(RQ)'] = row_medians

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
            # 只有一个 batch，跳过
            df.loc[method, ct] = np.nan
            continue

        # 取 embedding
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

        # 如果样本数太少（例如 <=1），直接跳过
        if n_samples <= 31:
            df.loc[method, ct] = np.nan
            continue

        # 先直接尝试计算 SAS
        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[method, ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            # 如果报错与 n_neighbors 有关，尝试用较小的邻居数重试（若函数支持该参数）
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                # 选择一个安全的 k：至少 1，最多 30（或 n_samples-1）
                k_try = max(1, min(30, n_samples - 1))
                try:
                    # 仅在该函数接受 n_neighbors 参数时有效
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[method, ct] = sas_tmp
                    print(f"{method}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{method}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[method, ct] = np.nan
            else:
                # 其它错误直接记录 NaN 并打印原因
                print(f"Error computing SAS for method {method}, celltype {ct}: {e}")
                df.loc[method, ct] = np.nan

row_medians = df.median(axis=1, skipna=True)
metric_df['CSAS(RQ)'] = row_medians

r_df = pyreadr.read_r(out_dir + "lisi_score_modality.rds")[None]
if metric_df.index.equals(r_df.index):
    metric_df['CiLISI(RQ)'] = r_df['CiLISI']
else:
    r_df2 = r_df.reindex(metric_df.index)
    metric_df['CiLISI(RQ)'] = r_df2['CiLISI']

metric_df['Bio conservation'] = (metric_df['KMeans NMI']+metric_df['KMeans ARI']+metric_df['Silhouette label']+metric_df['cLISI'])/4
metric_df['Batch correction'] = (metric_df['Silhouette batch']+metric_df['iLISI']+metric_df['KBET']+metric_df['Graph connectivity'])/4
metric_df['Total'] = metric_df['Bio conservation']*0.6 + metric_df['Batch correction']*0.4
metric_df['Modality mixing'] = (metric_df['CkBET']+metric_df['CSAS']+metric_df['CiLISI'])/3  
metric_df['RQ mixing'] = (metric_df['CkBET(RQ)']+metric_df['CSAS(RQ)']+metric_df['CiLISI(RQ)'])/3  

metric_df.loc['Metric Type'] = ['Bio conservation','Bio conservation','Bio conservation','Bio conservation',
                                'Batch correction','Batch correction','Batch correction','Batch correction',
                                'Modality mixing','Modality mixing','Modality mixing',
                                'RQ mixing','RQ mixing','RQ mixing',
                               'Aggregate score','Aggregate score','Aggregate score',
                               'Aggregate score','Aggregate score']

latent_data = metric_df
# 如果 latent_data 是一个 NumPy 数组，将其转换为 DataFrame
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    # 如果 latent_data 已经是 DataFrame，则直接使用
    latent_df = latent_data

# 保存为 CSV 文件
latent_df.to_csv(out_dir + 'metrics_new.csv', index=True)