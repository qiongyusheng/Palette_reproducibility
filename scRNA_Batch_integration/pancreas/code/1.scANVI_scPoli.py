import scanpy as sc
import scvi
import pyreadr
import anndata as ad
import numpy as np
import pandas as pd
from scvi.model.utils import mde
from sklearn.metrics import classification_report
from scarches.models.scpoli import scPoli
import warnings
warnings.filterwarnings('ignore')

adata = sc.read('./adata.h5ad')

early_stopping_kwargs = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}

data_dir = "human_pancreas"
condition_key = 'Batch'
cell_type_key = 'CellType'
npcs = 20

scvi.model.SCVI.setup_anndata(adata, batch_key=condition_key)
vae = scvi.model.SCVI(adata, n_layers=3, n_latent=20, gene_likelihood="nb")
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    labels_key=cell_type_key,
    unlabeled_category="Unknown",
)
lvae.train(max_epochs=20, n_samples_per_label=100)
res_scanvi = pd.DataFrame(lvae.get_latent_representation(adata))
res_scanvi.to_csv("./scanvi.csv",index=False)

adata.X = adata.X.astype(np.float32)
scpoli_model = scPoli(
    adata=adata,
    condition_keys=condition_key,
    cell_type_keys=cell_type_key,
    latent_dim=20,
    embedding_dims=5,
    recon_loss='nb',
)

scpoli_model.train(
    n_epochs=100,
    pretraining_epochs=95,
    early_stopping_kwargs=early_stopping_kwargs,
    eta=5,
)
res_scpoli = pd.DataFrame(scpoli_model.get_latent(adata))
res_scpoli.to_csv("./scpoli.csv", index=False)