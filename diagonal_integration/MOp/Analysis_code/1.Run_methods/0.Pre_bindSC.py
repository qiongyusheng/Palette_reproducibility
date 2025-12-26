import numpy as np
import pyreadr
import scipy.sparse as sp
from scipy.sparse.linalg import svds
import anndata as ad
import pandas as pd
import scanpy as sc

atac_co = sp.csr_matrix(pyreadr.read_r('./MERFISH_ATAC/ATAC/co.rds')[None])
adata = ad.AnnData(
    atac_co.T
)
del atac_co
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

rna = sp.csr_matrix(pyreadr.read_r('./MERFISH_ATAC/MERFISH/co.rds')[None])
matrix_product = adata.X @ rna
del adata
del rna
U, S, Vt = svds(matrix_product, k = 30)
del matrix_product
del S
# 保存转置后的U矩阵为CSV
np.savetxt('./MERFISH_ATAC/ATAC/bdsc30.csv', U.T, delimiter=',')
# 保存Vt矩阵为CSV
np.savetxt('./MERFISH_ATAC/MERFISH/bdsc30.csv', Vt, delimiter=',')
del U
del Vt

# 200
atac_co = sp.csr_matrix(pyreadr.read_r('./MERFISH_ATAC/ATAC/co200_RF.rds')[None])
adata = ad.AnnData(
    atac_co.T
)
del atac_co
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

rna = sp.csr_matrix(pyreadr.read_r('./MERFISH_ATAC/MERFISH/co200_RF.rds')[None])
matrix_product = adata.X @ rna
del adata
del rna
U, S, Vt = svds(matrix_product, k = 30)
del matrix_product
del S

# 保存转置后的U矩阵为CSV
np.savetxt('./MERFISH_ATAC/ATAC/bdsc30_200RF.csv', U.T, delimiter=',')
# 保存Vt矩阵为CSV
np.savetxt('./MERFISH_ATAC/MERFISH/bdsc30_200RF.csv', Vt, delimiter=',')
del U
del Vt

# 150
atac_co = sp.csr_matrix(pyreadr.read_r('./MERFISH_ATAC/ATAC/co150_RF.rds')[None])
adata = ad.AnnData(
    atac_co.T
)
del atac_co
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

rna = sp.csr_matrix(pyreadr.read_r('./MERFISH_ATAC/MERFISH/co150_RF.rds')[None])

matrix_product = adata.X @ rna
del adata
del rna
U, S, Vt = svds(matrix_product, k = 30)
del matrix_product
del S

# 保存转置后的U矩阵为CSV
np.savetxt('./MERFISH_ATAC/ATAC/bdsc30_150RF.csv', U.T, delimiter=',')
# 保存Vt矩阵为CSV
np.savetxt('./MERFISH_ATAC/MERFISH/bdsc30_150RF.csv', Vt, delimiter=',')
del U
del Vt

# 100
atac_co = sp.csr_matrix(pyreadr.read_r('./MERFISH_ATAC/ATAC/co100_RF.rds')[None])
adata = ad.AnnData(
    atac_co.T
)
del atac_co
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

rna = sp.csr_matrix(pyreadr.read_r('./MERFISH_ATAC/MERFISH/co100_RF.rds')[None])

matrix_product = adata.X @ rna
del adata
del rna
U, S, Vt = svds(matrix_product, k = 30)
del matrix_product
del S

# 保存转置后的U矩阵为CSV
np.savetxt('./MERFISH_ATAC/ATAC/bdsc30_100RF.csv', U.T, delimiter=',')
# 保存Vt矩阵为CSV
np.savetxt('./MERFISH_ATAC/MERFISH/bdsc30_100RF.csv', Vt, delimiter=',')
del U
del Vt

# 50
atac_co = sp.csr_matrix(pyreadr.read_r('./MERFISH_ATAC/ATAC/co50_RF.rds')[None])
adata = ad.AnnData(
    atac_co.T
)
del atac_co
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

rna = sp.csr_matrix(pyreadr.read_r('./MERFISH_ATAC/MERFISH/co50_RF.rds')[None])

matrix_product = adata.X @ rna
del adata
del rna
U, S, Vt = svds(matrix_product, k = 30)
del matrix_product
del S

# 保存转置后的U矩阵为CSV
np.savetxt('./MERFISH_ATAC/ATAC/bdsc30_50RF.csv', U.T, delimiter=',')
# 保存Vt矩阵为CSV
np.savetxt('./MERFISH_ATAC/MERFISH/bdsc30_50RF.csv', Vt, delimiter=',')
del U
del Vt