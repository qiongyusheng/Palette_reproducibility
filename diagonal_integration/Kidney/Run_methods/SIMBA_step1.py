import numpy as np
import pandas as pd
from scipy.io import mmread
import scanpy as sc
import anndata as ad
import json
import os
import pyreadr
import os
import simba as si
si.__version__

adata_CP = sc.read('/home/server/sqy/data/mosaic/RAW/Human_kidney/Muto-2021-ATAC.h5ad')
adata_CG = sc.read('/home/server/sqy/data/mosaic/RAW/Human_kidney/Muto-2021-RNA.h5ad')

adata_CP.obs_names = adata_CP.obs_names + '_ATAC'
adata_CG.obs_names = adata_CG.obs_names + '_RNA'

rds_rna = pyreadr.read_r('./Human_kidney/RNA/X.rds')
feat = rds_rna[None].index
feat = list(feat.astype(str)) if hasattr(feat, "astype") else list(map(str, feat))
gene_ids = pd.Index(adata_CG.var['gene_ids'].astype(str))   # 关键：转成 Index
idx = gene_ids.get_indexer(feat)   # Index 方法
present_mask = idx >= 0
idx_keep = idx[present_mask]
adata_CG = adata_CG[:, idx_keep].copy()
del rds_rna
del feat
del gene_ids
del present_mask
del idx
del idx_keep

si.pp.filter_peaks(adata_CP,min_n_cells=0)
si.pp.cal_qc_atac(adata_CP)
si.pl.violin(adata_CP,list_obs=['n_counts','n_peaks','pct_peaks'], list_var=['n_cells'])
si.pl.hist(adata_CP,list_obs=['n_counts','n_peaks','pct_peaks'], log=True, list_var=['n_cells'])
si.pp.pca(adata_CP, n_components=50)

si.pp.filter_genes(adata_CG,min_n_cells=0)
si.pp.cal_qc_rna(adata_CG)
si.pp.normalize(adata_CG,method='lib_size')
si.pp.log_transform(adata_CG)
si.pp.select_variable_genes(adata_CG, n_top_genes=2000)
si.pl.variable_genes(adata_CG,show_texts=True)

si.tl.discretize(adata_CG,n_bins=5)
si.pl.discretize(adata_CG,kde=False)
adata_CP.var['chr'] = adata_CP.var['chrom']
adata_CP.var['start'] = adata_CP.var['chromStart']
adata_CP.var['end'] = adata_CP.var['chromEnd']
adata_CG_atac = si.tl.gene_scores(adata_CP,genome='hg38',use_gene_weigt=True, use_top_pcs=False)
si.pp.filter_genes(adata_CG_atac,min_n_cells=3)
si.pp.cal_qc_rna(adata_CG_atac)
si.pp.normalize(adata_CG_atac,method='lib_size')
si.pp.log_transform(adata_CG_atac)

adata_CrnaCatac = si.tl.infer_edges(adata_CG, adata_CG_atac, n_components=15, k=15)

si.pl.node_similarity(adata_CrnaCatac,cutoff=0.5)
si.tl.trim_edges(adata_CrnaCatac, cutoff=0.5)

si.tl.gen_graph(list_CP=[adata_CP],
                list_CG=[adata_CG],
                list_CC=[adata_CrnaCatac],
                copy=False,
                use_highly_variable=True,
                use_top_pcs=False,
                dirname='graph0')
dict_config = si.settings.pbg_params.copy()
dict_config['dimension'] = 20
si.tl.pbg_train(pbg_params = dict_config, auto_wd=True, save_wd=True, output='model')
dict_adata = si.read_embedding()
adata_all = si.tl.embed(adata_ref=dict_adata['C2'],
                        list_adata_query=[dict_adata['C']],
                        use_precomputed=False)
latent_data = adata_all.X

if isinstance(latent_data, np.ndarray):
    latent_df = pd.DataFrame(latent_data)
else:
    latent_df = latent_data

latent_df.to_csv('/home/server/sqy/data/mosaic/result/NC_revise_1/weak_link/kidney/simba.csv', index=False)

adata_all.obs_names.to_series().to_csv("/home/server/sqy/data/mosaic/result/NC_revise_1/weak_link/kidney/simba_cell_name.csv", 
                                       index=False, header=False)