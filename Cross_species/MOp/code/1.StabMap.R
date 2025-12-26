library(StabMap)
library(scran)
library(Seurat)
library(patchwork)

setwd('./cross_species')

# human
# 1.SNARE (RNA+ATAC)
h_s_atac <- as(as.matrix(readRDS('./human/SNARE_Lein/atac.rds')),'dgCMatrix')
h_s_rna <- as(as.matrix(readRDS('./human/SNARE_Lein/co.rds')),'dgCMatrix')
h_s_m <- readRDS('./human/SNARE_Lein/m.rds')

# 2.10XV3 RNA
h_10xv3_rna <- as(as.matrix(readRDS('./human/10XV3_Lein/co.rds')),'dgCMatrix')
h_10xv3_m <- readRDS('./human/10XV3_Lein/m.rds')

# 3. 10XMultiome ATAC
h_10xmulti_atac <- as(as.matrix(readRDS('./human/Multiome_Zemke/atac.rds')),'dgCMatrix')
h_10xmulti_m <- readRDS('./human/Multiome_Zemke/m.rds')

rownames(h_10xmulti_atac) <- rownames(h_s_atac) <- paste0('human_',rownames(h_s_atac))

# mouse
# 1.10XMultiome (RNA+ATAC)
m_10xmulti_rna <- as(as.matrix(readRDS('./mouse/10XMultiome_Zemke/co.rds')),'dgCMatrix')
m_10xmulti_atac <- as(as.matrix(readRDS('./mouse/10XMultiome_Zemke/atac.rds')),'dgCMatrix')
m_10xmulti_m <- readRDS('./mouse/10XMultiome_Zemke/m.rds')

# 2. 10XV3 RNA
m_10xv3_rna <- as(as.matrix(readRDS('./mouse/10XV3_Lein/co.rds')),'dgCMatrix')
m_10xv3_m <- readRDS('./mouse/10XV3_Lein/m.rds')

# 3. 10X ATAC
# m_10x_atac <- as.matrix(readRDS('./mouse/10X/SVD.rds'))
m_10x_atac <- as(as.matrix(readRDS('./mouse/10X/atac.rds')),'dgCMatrix')
m_10x_m <- readRDS('./mouse/10X/m.rds')

# Marmoset
# 1.SNARE (RNA+ATAC)
mar_snare_rna <- as(as.matrix(readRDS('./Marmoset/SNARE_Lein/co.rds')),'dgCMatrix')
mar_snare_m <- readRDS('./Marmoset/SNARE_Lein/m.rds')

# 2.10XV3 RNA
mar_10xv3_rna <- as(as.matrix(readRDS('./Marmoset/10XV3_Lein/co.rds')),'dgCMatrix')
mar_10xv3_m <- readRDS('./Marmoset/10XV3_Lein/m.rds')

h_s_atac <-  NormalizeData(h_s_atac)
h_s_rna <-  NormalizeData(h_s_rna)

h_10xv3_rna <-  NormalizeData(h_10xv3_rna)

h_10xmulti_atac <-  NormalizeData(h_10xmulti_atac)

m_10xmulti_rna <- NormalizeData(m_10xmulti_rna)
m_10xmulti_atac <- NormalizeData(m_10xmulti_atac)

m_10xv3_rna <- NormalizeData(m_10xv3_rna)

m_10x_atac <- NormalizeData(m_10x_atac)

mar_snare_rna <- NormalizeData(mar_snare_rna)

mar_10xv3_rna <- NormalizeData(mar_10xv3_rna)

assay_list = list(human1 = rbind(h_s_atac,h_s_rna),
                  human2 = h_10xv3_rna,
                  human3 = h_10xmulti_atac,
                  
                  mouse1 = rbind(m_10xmulti_rna,m_10xmulti_atac),
                  mouse2 = m_10xv3_rna,
                  mouse3 = m_10x_atac,
                 
                  marmoset1 = mar_snare_rna,
                  marmoset2 = mar_10xv3_rna)

rm(h_s_atac,h_s_rna,h_10xv3_rna,h_10xmulti_atac,
   m_10xmulti_rna,m_10xmulti_atac,m_10xv3_rna,m_10x_atac,
  mar_snare_rna,mar_10xv3_rna)

mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("human1"),
               maxFeatures = 1000000,
              ncomponentsReference = 20,
              ncomponentsSubset = 20)

meta <- readRDS('./Cross_species/MOp/output/meta.rds')

stab <- stab[rownames(meta),]
out_dir <- './Cross_species/MOp/output'
write.csv(as.matrix(stab),
          paste0(out_dir,'/stab.csv'),
          row.names = FALSE)
