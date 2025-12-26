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


dataset = sc.read('./adata.h5ad')

meta = pyreadr.read_r("./meta.rds")[None]

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
cell_type_key = ['label_miss5s','label_miss10s','label_miss20s']
mod = ['_miss5','_miss10','_miss20']
npcs = 20

for i in [0,1,2]:
    # ----------- Data -------------
    adata = ad.AnnData(
      X=dataset.X,
      obs=meta.copy()
    )
    adata.obs_names = dataset.obs_names
    adata.var_names = dataset.var_names
    adata.X = adata.X.astype(np.float32)

    # ----------- SCANVI -------------
    scvi.model.SCVI.setup_anndata(adata, batch_key=condition_key)
    vae = scvi.model.SCVI(
        adata,
        n_layers=3,
        n_latent=20,
        gene_likelihood="nb"
    )
    lvae = scvi.model.SCANVI.from_scvi_model(
        vae,
        adata=adata,
        labels_key=cell_type_key[i],
        unlabeled_category="Unknown",
    )
    lvae.train(
        max_epochs=20,
        n_samples_per_label=100,
    )

    res_scanvi = pd.DataFrame(lvae.get_latent_representation(adata))
    res_scanvi.to_csv(
        f"./label_miss/scanvi{mod[i]}.csv",
        index=False
    )

    # ----------- scPoli -------------
    scpoli_model = scPoli(
        adata=adata,
        condition_keys=condition_key,
        cell_type_keys=cell_type_key[i],
        unknown_ct_names = ["Unknown"],
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
    res_scpoli.to_csv(
        f"./label_miss/scpoli{mod[i]}.csv",
        index=False
    )
