import scanpy as sc
import scvi
import pyreadr
import anndata
import numpy as np
import pandas as pd
from scvi.model.utils import mde
from sklearn.metrics import classification_report
from scarches.models.scpoli import scPoli
import warnings
warnings.filterwarnings('ignore')

data_dir = "WAT"
condition_key = 'Batch'
cell_type_key = 'CellType'
npcs = 20
early_stopping_kwargs = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}
########################################## Run on HVGs
# prepare adata object
dataset1 = pyreadr.read_r("./human_mouse_WAT/human/rna_homo.rds")[None]
dataset2 = pyreadr.read_r("./human_mouse_WAT/mouse/rna_homo.rds")[None]
dataset = pd.concat([dataset1, dataset2], axis=1)
del dataset1
del dataset2
meta = pyreadr.read_r("./meta.rds")[None]
adata = anndata.AnnData(X=dataset.T, obs=meta)
adata.X = adata.X.astype(np.float32)

adata = anndata.AnnData(X=dataset.T, obs=meta)
adata.X = adata.X.astype(np.float32)
scvi.model.SCVI.setup_anndata(adata, batch_key=condition_key)
vae = scvi.model.SCVI(adata, n_layers=3, n_latent=30, gene_likelihood="nb")
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    labels_key=cell_type_key,
    unlabeled_category="Unknown",
)
lvae.train(max_epochs=20, n_samples_per_label=100)#set batch size for this dataset , batch_size = 100
res_scanvi = pd.DataFrame(lvae.get_latent_representation(adata))
res_scanvi.to_csv("./scanvi.csv",index=False)

scpoli_model = scPoli(
    adata=adata,
    condition_keys=condition_key,
    cell_type_keys=cell_type_key,
    latent_dim=25,
    embedding_dims=5,
    recon_loss='nb',
)

scpoli_model.train(
    n_epochs=100,
    pretraining_epochs=95,
    early_stopping_kwargs=early_stopping_kwargs,
    eta=5,
)

latent = scpoli_model.get_latent(adata)
if isinstance(latent, np.ndarray):
    latent = latent[:, :20]
    latent_df = pd.DataFrame(latent)
else:
    latent_df = latent.iloc[:, :20]

latent_df.to_csv('./scpoli.csv', index=False)