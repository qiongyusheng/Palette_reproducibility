import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.io as sio
from scipy.sparse import csr_matrix, issparse
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.decomposition import TruncatedSVD
from sklearn import metrics
import muon
import episcanpy as epi

# ---------- helpers ----------

def cluster_lov(data, n):
    n = int(n)
    epi.tl.getNClusters(data, n, method='louvain')
    return data.obs['louvain'].astype(str).to_numpy()

def cluster_len(data, n):
    n = int(n)
    epi.tl.getNClusters(data, n, method='leiden')
    return data.obs['leiden'].astype(str).to_numpy()

def run_lsi(X, n_components=50):
    """
    X: csr_matrix / ndarray (cells × peaks)
    return: ndarray (cells × n_components)
    """
    # SVD
    svd = TruncatedSVD(n_components=n_components, random_state=42)
    Z = svd.fit_transform(X)
    return Z

def ensure_2d_csr(x):
    if issparse(x):
        return x.tocsr()
    arr = np.asarray(x)
    if arr.ndim != 2:
        raise ValueError("Input must be 2D.")
    return csr_matrix(arr)

import os
import pandas as pd
import scipy.io as sio
from scipy.sparse import csr_matrix

def read_matrix(method_name: str, file_path: str, is_atac: bool, csv_header='infer'):
    """
    返回 cells × features 的 csr_matrix，读取/转置逻辑按你的原始代码：
    - is_atac = True:
        Palette        -> 读取 .mtx 并转置
        MIDAS, scVAEIT -> 读取 .csv（不转置）
        其他            -> 读取 .mtx（不转置）
    - is_atac = False:
        Palette                     -> 读取 .csv 并转置
        MIDAS, Multigrate, scVAEIT  -> 读取 .csv（不转置）
    其余情况：默认按 .mtx（不转置）或 .csv（不转置）读取
    """
    fname = os.path.basename(file_path)
    is_mtx = fname.endswith(".mtx") or fname.endswith(".mtx.gz")
    method_name = str(method_name)

    def _read_csv(path, transpose=False):
        # header='infer' 会在首行非数值时自动当作表头，若你明确无表头可传 header=None
        df = pd.read_csv(path, header=csv_header)
        arr = df.to_numpy()
        if transpose:
            arr = arr.T
        return csr_matrix(arr)

    def _read_mtx(path, transpose=False):
        M = sio.mmread(path).tocsr()
        return M.T if transpose else M

    if is_atac:
        # ----- ATAC 分支：按你原始逻辑 -----
        if method_name == 'Palette':
            if not is_mtx:
                raise ValueError(f"[ATAC/Palette] 期望 MTX 文件，但得到: {file_path}")
            X = _read_mtx(file_path, transpose=True)          # 转置
        elif method_name in ('MIDAS', 'scVAEIT'):
            if is_mtx:
                raise ValueError(f"[ATAC/{method_name}] 期望 CSV 文件，但得到 MTX: {file_path}")
            X = _read_csv(file_path, transpose=False)         # 不转置
        else:
            # 其他方法 -> 读 mtx，不转置
            if not is_mtx:
                # 若你确实给了 CSV，这里也可兼容地读 CSV
                X = _read_csv(file_path, transpose=False)
            else:
                X = _read_mtx(file_path, transpose=False)
    else:
        # ----- RNA 分支：按你原始逻辑 -----
        if method_name == 'Palette':
            if is_mtx:
                raise ValueError(f"[RNA/Palette] 期望 CSV 文件，但得到 MTX: {file_path}")
            X = _read_csv(file_path, transpose=True)          # 转置
        elif method_name in ('MIDAS', 'Multigrate', 'scVAEIT'):
            if is_mtx:
                raise ValueError(f"[RNA/{method_name}] 期望 CSV 文件，但得到 MTX: {file_path}")
            X = _read_csv(file_path, transpose=False)         # 不转置
        else:
            # 其他方法：默认优先按 MTX 读
            X = _read_mtx(file_path, transpose=False) if is_mtx else _read_csv(file_path, transpose=False)

    # 最终确保是 cells × features
    if X.shape[0] <= 1 or X.shape[1] <= 1:
        # 简单 sanity check；必要时你可在这里加更多断言
        pass
    return X


# ----------- main ------------

def cal_metrics(
    path,
    data_path,              # list[str]：每个方法的数据文件相对 path 的路径
    method_name,            # list[str]：方法名（与 data_path 对齐）
    meta_path,              # RDS/CSV：包含真实标签
    label_key,              # 真实标签列名
    norm,                   # list[str|None]：每个方法的标准化方法
    dim_method,             # 'pca' | 'lsi'
    dim,                    # int：降维维度
    is_ATAC=False,
    clust_method='louvain'
):
    assert len(method_name) == len(data_path) == len(norm), "参数长度必须一致"

    # 真实标签
    if meta_path.endswith(".rds"):
        import pyreadr
        meta = pyreadr.read_r(meta_path)[None]
    else:
        meta = pd.read_csv(meta_path)
    label = meta[label_key].astype(str).to_numpy()

    metric_df = pd.DataFrame(
        index=method_name,
        columns=['ARI','AMI','FMI','ASW','homogeneity score','completeness score','NMI'],
        dtype=float
    )

    for i, mname in enumerate(method_name):
        f = os.path.join(path, data_path[i])

        # -------- 读取矩阵（拿到 cells × features）--------
        # 你原始代码中对不同方法不同转置很混乱，这里给一个统一入口：
        # 如果你已确定某些文件一定是 features × cells，可以给这些文件单独加个标志位；
        # 这里简单做法：先读，再按需要转置到 cells × features。
        # X = read_matrix(f)        # csr_matrix(cells × features)
        X = read_matrix(method_name[i], os.path.join(path, data_path[i]), is_ATAC, csv_header='infer')
        n_cells = X.shape[0]

        # 与标签长度检查（强烈建议一致）
        if len(label) != n_cells:
            raise ValueError(f"[{mname}] 标签长度 {len(label)} 与矩阵行数 {n_cells} 不一致！")

        adata = ad.AnnData(X)

        # -------- 标准化 --------
        if norm[i] == 'log':
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
        elif norm[i] == 'clr':
            # 通常 CLR 用在 ADT；ATAC 用 TF-IDF；RNA 用 normalize+log1p
            muon.prot.pp.clr(adata)
        elif norm[i] == 'tfidf':
            tfidf = TfidfTransformer(norm='l2', use_idf=True, smooth_idf=True, sublinear_tf=False)
            adata.X = tfidf.fit_transform(adata.X)

        # -------- 降维 --------
        if dim_method == 'pca':
            # 对稀疏矩阵先转为密集会很耗内存；小心使用
            sc.pp.scale(adata)          # 需要 dense；若非常大，请改用 IncrementalPCA
            sc.tl.pca(adata, n_comps=dim)
            Z = adata.obsm['X_pca']     # (n_cells, dim)
        elif dim_method == 'lsi':
            Z = run_lsi(adata.X, n_components=dim)
            adata.obsm['X_pca'] = Z     # 统一用 'X_pca' 作为后续度量入口
        else:
            raise ValueError("dim_method must be 'pca' or 'lsi'")

        # -------- 聚类 --------
        n_clusters = len(np.unique(label))
        sc.pp.neighbors(adata)
        if clust_method == 'louvain':
            pred = cluster_lov(adata, n_clusters)
        elif clust_method == 'leiden':
            pred = cluster_len(adata, n_clusters)
        else:
            raise ValueError("clust_method must be 'louvain' or 'leiden'")

        # -------- 指标 --------
        metric_df.loc[mname, 'ARI'] = metrics.adjusted_rand_score(label, pred)
        metric_df.loc[mname, 'AMI'] = metrics.adjusted_mutual_info_score(label, pred)
        metric_df.loc[mname, 'FMI'] = metrics.fowlkes_mallows_score(label, pred)

        # silhouette：用嵌入空间欧氏距离（推荐）；你之前的 1-corr 也可，但要小心 NaN
        try:
            asw = metrics.silhouette_score(Z, label, metric='euclidean')
        except Exception:
            # 退而求其次，用相关距离（可能出现 NaN）
            C = np.corrcoef(Z.T) if Z.shape[0] < Z.shape[1] else np.corrcoef(Z)
            D = 1 - C
            asw = metrics.silhouette_score(D, label, metric='precomputed')
        metric_df.loc[mname, 'ASW'] = asw

        metric_df.loc[mname, 'homogeneity score']  = metrics.homogeneity_score(label, pred)
        metric_df.loc[mname, 'completeness score'] = metrics.completeness_score(label, pred)
        metric_df.loc[mname, 'NMI']                = metrics.normalized_mutual_info_score(label, pred)

    metric_df['Method'] = method_name
    return metric_df
