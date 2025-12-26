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

data_dir = './Benchmarking/Supervised_using_modified_label/Label_miss/output/BMMC_s2/'

data_base = ['Palette_34','Palette_35','Palette_36','Palette_24','Palette_26']
meta_base = ['meta_34','meta_35','meta_36','meta_24','meta_26']
mod = ['_miss5','_miss10','_miss20']
data_name = [d + m for d in data_base for m in mod]
meta_name = [m for m in meta_base for _ in range(3)]

out_dir = './Benchmarking/Supervised_using_modified_label/Label_miss/output/BMMC_s2/'
meta_dir = data_dir

metric_df = pd.DataFrame(np.nan, index=data_name, columns=['KMeans NMI','KMeans ARI',
                                                           'Silhouette label','cLISI',
                                                           'Silhouette batch','iLISI',
                                                           'KBET','Graph connectivity',
                                                           'CkBET','CSAS','CiLISI'])

for i in range(len(meta_name)):
    meta = pyreadr.read_r(meta_dir + meta_name[i]+ ".rds")
    meta = meta[None]
    dataset = np.random.rand(meta.shape[0], 2)
    adata = anndata.AnnData(X=dataset, obs=meta)
    res = pd.read_csv(data_dir + data_name[i] + ".csv")
    adata.obsm[data_name[i]] = res.values

    condition_key = 'batch'
    cell_type_key = 'celltype.l2'
    batchcorr = BatchCorrection(pcr_comparison = False)
    bioconser = BioConservation(isolated_labels=False)
    bm = Benchmarker(
        adata,
        batch_key=condition_key,
        label_key=cell_type_key,
        embedding_obsm_keys=list([data_name[i]]),
        # pre_integrated_embedding_obsm_key = 'Raw',
        bio_conservation_metrics = bioconser,
        batch_correction_metrics = batchcorr,
        n_jobs=6,
        )
    bm.benchmark()
    df = bm.get_results(min_max_scale=False)
    metric_df.loc[data_name[i],metric_df.columns[0:8]] = df.loc[data_name[i],df.columns[0:8]]

    condition_key = 'modality'
    cell_type_key = 'celltype.l2'

    if pd.api.types.is_categorical_dtype(adata.obs[cell_type_key]):
        ct_list = list(adata.obs[cell_type_key].cat.categories)
    else:
        ct_list = sorted(adata.obs[cell_type_key].unique().tolist())
    # 4) 初始化结果表（method x celltype）
    df = pd.DataFrame(np.nan, index=list([data_name[i]]), columns=ct_list)

    batchcorr = BatchCorrection(pcr_comparison = False,
                            graph_connectivity=False,
                            ilisi_knn=False,
                            silhouette_batch=False)
    bioconser = BioConservation(isolated_labels = False,
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
        embedding_obsm_keys=list([data_name[i]]),
        bio_conservation_metrics=bioconser,
        batch_correction_metrics=batchcorr,
        n_jobs=6,)
        bm.benchmark()
        tmp = bm.get_results(min_max_scale=False)
         # 兼容性赋值（若 tmp 中缺少 method 或 KBET 列则填 NaN）
        for method in list([data_name[i]]):
            if method in tmp.index and 'KBET' in tmp.columns:
                df.loc[method, ct] = tmp.loc[method, 'KBET']
            else:
                df.loc[method, ct] = np.nan

    row_medians = df.median(axis=1, skipna=True)
    metric_df.at[data_name[i], 'CkBET'] = row_medians

    df = pd.DataFrame(np.nan, index=list([data_name[i]]), columns=ct_list)

    for ct in ct_list:
        adata_ct = adata[adata.obs[cell_type_key] == ct].copy()
        if adata_ct.n_obs == 0:
            df.loc[data_name[i], ct] = np.nan
            continue

        unique_batches = adata_ct.obs[condition_key].unique()
        if len(unique_batches) < 2:
            # 只有一个 batch，跳过
            df.loc[data_name[i], ct] = np.nan
            continue

        # 取 embedding
        emb = adata_ct.obsm[data_name[i]]
        if emb.shape[0] != adata_ct.n_obs:
            if emb.shape[1] == adata_ct.n_obs:
                emb = emb.T
            else:
                raise ValueError(
                    f"Embedding shape {emb.shape} 不匹配 adata_ct.n_obs={adata_ct.n_obs} "
                    f"（method={[data_name[i]]}, celltype={ct}）。"
                )

        modal = np.array(adata_ct.obs[condition_key])
        n_samples = adata_ct.n_obs

        # 如果样本数太少（例如 <=1），直接跳过
        if n_samples <= 31:
            df.loc[data_name[i], ct] = np.nan
            continue

        # 先直接尝试计算 SAS
        try:
            sas_tmp = seurat_alignment_score(emb, modal)
            df.loc[data_name[i], ct] = sas_tmp
        except Exception as e:
            msg = str(e)
            # 如果报错与 n_neighbors 有关，尝试用较小的邻居数重试（若函数支持该参数）
            if 'n_neighbors' in msg or 'n_samples' in msg or 'Expected n_neighbors' in msg:
                # 选择一个安全的 k：至少 1，最多 30（或 n_samples-1）
                k_try = max(1, min(30, n_samples - 1))
                try:
                    # 仅在该函数接受 n_neighbors 参数时有效
                    sas_tmp = seurat_alignment_score(emb, modal, n_neighbors=k_try)
                    df.loc[[data_name[i]], ct] = sas_tmp
                    print(f"{[data_name[i]]}, {ct}: retried with n_neighbors={k_try} succeeded")
                except Exception as e2:
                    print(f"{[data_name[i]]}, {ct}: retry with n_neighbors={k_try} failed: {e2}")
                    df.loc[[data_name[i]], ct] = np.nan
            else:
                # 其它错误直接记录 NaN 并打印原因
                print(f"Error computing SAS for method {[data_name[i]]}, celltype {ct}: {e}")
                df.loc[[data_name[i]], ct] = np.nan
    row_medians = df.median(axis=1, skipna=True)
    metric_df.at[data_name[i], 'CSAS'] = row_medians

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
# 如果 latent_data 是一个 NumPy 数组，将其转换为 DataFrame
if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    # 如果 latent_data 已经是 DataFrame，则直接使用
    latent_df = latent_data

# 保存为 CSV 文件
latent_df.to_csv(out_dir + 'metrics.csv', index=True)