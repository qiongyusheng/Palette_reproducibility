
# Downlaod Link: https://explore.data.humancellatlas.org/projects/cc95ff89-2e68-4a08-a234-480eca21ce79/project-matrices

import loompy
import scipy.sparse as sp
import numpy as np
import os
from scipy.io import mmwrite


loom_path = "./1M-immune-human-immune-10XV2.loom"
ds = loompy.connect(loom_path, mode="r")
blocks = []
chunk = 2000

for i in range(0, ds.shape[0], chunk):
    print(f"Reading rows {i} - {min(i+chunk, ds.shape[0])}")
    X_block = ds[i:i+chunk, :]          # small dense block (genes x cells)
    blocks.append(sp.csr_matrix(X_block))

X = sp.vstack(blocks).tocsr()
print("Sparse matrix shape:", X.shape)

gene_key = "Gene" if "Gene" in ds.ra.keys() else list(ds.ra.keys())[0]
barcode_key = "CellID" if "CellID" in ds.ca.keys() else list(ds.ca.keys())[0]

genes = ds.ra[gene_key].astype(str)
barcodes = ds.ca[barcode_key].astype(str)

output_dir = "./immune_mtx/"
os.makedirs(output_dir, exist_ok=True)

X_to_write = X.astype(np.int32)

mmwrite(os.path.join(output_dir, "matrix.mtx"), X_to_write)

with open(os.path.join(output_dir, "genes.tsv"), "w") as f:
    for g in genes:
        f.write(g + "\n")

with open(os.path.join(output_dir, "barcodes.tsv"), "w") as f:
    for bc in barcodes:
        f.write(bc + "\n")
ds.close()

print("Matrix Market (.mtx) export completed!")
