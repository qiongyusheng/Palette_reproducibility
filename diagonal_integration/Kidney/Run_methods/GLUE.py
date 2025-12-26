import anndata as ad
import pandas as pd
import networkx as nx
import scanpy as sc
import scglue
from matplotlib import rcParams
import pyreadr
import numpy as np
from scipy.sparse import csr_matrix
atac = sc.read('./Muto-2021-ATAC.h5ad')
rna = sc.read('./Muto-2021-RNA.h5ad')
rds_rna = pyreadr.read_r('./X.rds')
feat = rds_rna[None].index
feat = list(feat.astype(str)) if hasattr(feat, "astype") else list(map(str, feat))

gene_ids = pd.Index(rna.var['gene_ids'].astype(str))   # 关键：转成 Index
idx = gene_ids.get_indexer(feat)   # Index 方法
present_mask = idx >= 0
idx_keep = idx[present_mask]
rna = rna[:, idx_keep].copy()

del rds_rna
del feat
del gene_ids
del present_mask
del idx
del idx_keep
rna.layers["counts"] = rna.X.copy()
rna.var.loc[:, ["chrom", "chromStart", "chromEnd"]].head()
atac.var_names[:5]

graph = scglue.genomics.rna_anchored_prior_graph(rna, atac)
graph.number_of_nodes(), graph.number_of_edges()

# Graph node covers all omic features
all(graph.has_node(gene) for gene in rna.var_names), \
all(graph.has_node(peak) for peak in atac.var_names)

# Edge attributes contain weights and signs
for _, e in zip(range(5), graph.edges):
    print(f"{e}: {graph.edges[e]}")
# Each node has a self-loop
all(graph.has_edge(gene, gene) for gene in rna.var_names), \
all(graph.has_edge(peak, peak) for peak in atac.var_names)

# Graph is symmetric
all(graph.has_edge(j, i) for i, j, _ in graph.edges)

import anndata
import itertools
import networkx as nx
import pandas as pd
import scanpy as sc
import scglue
import seaborn as sns
from matplotlib import rcParams

rna.layers["counts"] = rna.X.copy()
sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3")
sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
sc.pp.scale(rna)
sc.tl.pca(rna, n_comps=100, svd_solver="auto")

scglue.models.configure_dataset(
    rna, "NB", use_highly_variable=True,
    use_layer="counts", use_rep="X_pca"
)

scglue.data.lsi(atac, n_components=100, n_iter=15)

scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=True,
    use_rep="X_lsi"
)

graph = graph.subgraph(itertools.chain(
    rna.var.query("highly_variable").index,
    atac.var.query("highly_variable").index
))

model = scglue.models.SCGLUEModel(
    {"rna": rna, "atac": atac}, graph,
    latent_dim=20
)


model.compile()

glue = model.fit(
    {"rna": rna, "atac": atac},
    graph
)

rna.obsm["X_glue"] = model.encode_data("rna", rna)
atac.obsm["X_glue"] = model.encode_data("atac", atac)

import pandas as pd

rna.obs_names  = rna.obs_names.astype(str)
atac.obs_names = atac.obs_names.astype(str)
rna.obs_names  = rna.obs_names.map(lambda s: f"rna{s}")
atac.obs_names = atac.obs_names.map(lambda s: f"atac{s}")
combined = anndata.concat([rna, atac])

glue_array = combined.obsm['X_glue']
import numpy as np
glue_array = pd.DataFrame(glue_array)
glue_array.index = np.concatenate([rna.obs_names,atac.obs_names])

glue_array.to_csv('./glue.csv', index=True)
