
meta <- readRDS('./Benchmarking/Unsupervised/TEA_s2/output/meta.rds')
head(meta)

meta$label_shuffled5 <- meta$celltype
meta$label_shuffled10 <- meta$celltype
meta$label_shuffled20 <- meta$celltype

sample_size <- ceiling(0.05 * nrow(meta))
set.seed(123)
idx <- sample(1:nrow(meta), sample_size,replace = F)
meta$label_shuffled5[idx] <- sample(meta$label_shuffled5[idx])

sample_size <- ceiling(0.1 * nrow(meta))
set.seed(123)
idx <- sample(1:nrow(meta), sample_size,replace = F)
meta$label_shuffled10[idx] <- sample(meta$label_shuffled10[idx])

sample_size <- ceiling(0.2 * nrow(meta))
set.seed(123)
idx <- sample(1:nrow(meta), sample_size,replace = F)
meta$label_shuffled20[idx] <- sample(meta$label_shuffled20[idx])

saveRDS(meta,'./Benchmarking/Supervised_using_modified_label/Label_shuffled/output/TEA_s1/meta.rds')
saveRDS(meta,'./Benchmarking/Supervised_using_modified_label/Label_shuffled/output/TEA_s2/meta.rds')



# BMMC_s1

meta <- readRDS('./Benchmarking/Unsupervised/BMMC_s1/output/meta.rds')

head(meta)

meta$label_shuffled5 <- meta$celltype.l2
meta$label_shuffled10 <- meta$celltype.l2
meta$label_shuffled20 <- meta$celltype.l2

sample_size <- ceiling(0.05 * nrow(meta))
set.seed(123)
idx <- sample(1:nrow(meta), sample_size,replace = F)
meta$label_shuffled5[idx] <- sample(meta$label_shuffled5[idx])

sample_size <- ceiling(0.1 * nrow(meta))
set.seed(123)
idx <- sample(1:nrow(meta), sample_size,replace = F)
meta$label_shuffled10[idx] <- sample(meta$label_shuffled10[idx])

sample_size <- ceiling(0.2 * nrow(meta))
set.seed(123)
idx <- sample(1:nrow(meta), sample_size,replace = F)
meta$label_shuffled20[idx] <- sample(meta$label_shuffled20[idx])

meta_drop1 <- dplyr::filter(meta,batch != 'CITE1')
meta_drop2 <- dplyr::filter(meta,batch != 'CITE2')
meta_drop4 <- dplyr::filter(meta,batch != 'Multiome1')
meta_drop5 <- dplyr::filter(meta,batch != 'Multiome2')
meta_drop6 <- dplyr::filter(meta,batch != 'Multiome3')

saveRDS(meta_drop1,
        './Benchmarking/Supervised_using_modified_label/Label_shuffled/output/BMMC_s1/meta_drop1.rds')
saveRDS(meta_drop2,
        './Benchmarking/Supervised_using_modified_label/Label_shuffled/output/BMMC_s1/meta_drop2.rds')
saveRDS(meta_drop4,
        './Benchmarking/Supervised_using_modified_label/Label_shuffled/output/BMMC_s1/meta_drop4.rds')
saveRDS(meta_drop5,
        './Benchmarking/Supervised_using_modified_label/Label_shuffled/output/BMMC_s1/meta_drop5.rds')
saveRDS(meta_drop6,
        './Benchmarking/Supervised_using_modified_label/Label_shuffled/output/BMMC_s1/meta_drop6.rds')

meta_c1 <- dplyr::filter(meta,batch == 'CITE1')
meta_c2 <- dplyr::filter(meta,batch == 'CITE2')
meta_c3 <- dplyr::filter(meta,batch == 'CITE3')

meta_m1 <- dplyr::filter(meta,batch == 'Multiome1')
meta_m2 <- dplyr::filter(meta,batch == 'Multiome2')
meta_m3 <- dplyr::filter(meta,batch == 'Multiome3')

meta_BMMC_s <- readRDS('./Benchmarking/Unsupervised/BMMC_s2/output/meta_34.rds')

head(meta_BMMC_s)

unique(meta_BMMC_s$batch)

meta_BMMC_s$label_shuffled5 <- c(meta_c1$label_shuffled5,meta_c1$label_shuffled5,
                                meta_c2$label_shuffled5,meta_c2$label_shuffled5,
                                meta_c3$label_shuffled5,
                                meta_m1$label_shuffled5,
                                meta_m2$label_shuffled5,meta_m2$label_shuffled5,
                                meta_m3$label_shuffled5,meta_m3$label_shuffled5)

meta_BMMC_s$label_shuffled10 <- c(meta_c1$label_shuffled10,meta_c1$label_shuffled10,
                                meta_c2$label_shuffled10,meta_c2$label_shuffled10,
                                meta_c3$label_shuffled10,
                                meta_m1$label_shuffled10,
                                meta_m2$label_shuffled10,meta_m2$label_shuffled10,
                                meta_m3$label_shuffled10,meta_m3$label_shuffled10)

meta_BMMC_s$label_shuffled20 <- c(meta_c1$label_shuffled20,meta_c1$label_shuffled20,
                                meta_c2$label_shuffled20,meta_c2$label_shuffled20,
                                meta_c3$label_shuffled20,
                                meta_m1$label_shuffled20,
                                meta_m2$label_shuffled20,meta_m2$label_shuffled20,
                                meta_m3$label_shuffled20,meta_m3$label_shuffled20)

saveRDS(meta_BMMC_s,
        './Benchmarking/Supervised_using_modified_label/Label_shuffled/output/BMMC_s2/meta_34.rds')

meta_BMMC_s_24 <- readRDS('./Benchmarking/Unsupervised/BMMC_s2/output/meta_24.rds')

unique(meta_BMMC_s_24$batch)

meta_BMMC_s_24$label_shuffled5 <- c(meta_c1$label_shuffled5,meta_c1$label_shuffled5,
                                meta_c2$label_shuffled5,
                                meta_c3$label_shuffled5,meta_c3$label_shuffled5,
                                meta_m1$label_shuffled5,
                                meta_m2$label_shuffled5,meta_m2$label_shuffled5,
                                meta_m3$label_shuffled5,meta_m3$label_shuffled5)

meta_BMMC_s_24$label_shuffled10 <- c(meta_c1$label_shuffled10,meta_c1$label_shuffled10,
                                meta_c2$label_shuffled10,
                                meta_c3$label_shuffled10,meta_c3$label_shuffled10,
                                meta_m1$label_shuffled10,
                                meta_m2$label_shuffled10,meta_m2$label_shuffled10,
                                meta_m3$label_shuffled10,meta_m3$label_shuffled10)

meta_BMMC_s_24$label_shuffled20 <- c(meta_c1$label_shuffled20,meta_c1$label_shuffled20,
                                meta_c2$label_shuffled20,
                                meta_c3$label_shuffled20,meta_c3$label_shuffled20,
                                meta_m1$label_shuffled20,
                                meta_m2$label_shuffled20,meta_m2$label_shuffled20,
                                meta_m3$label_shuffled20,meta_m3$label_shuffled20)
saveRDS(meta_BMMC_s_24,
        './Benchmarking/Supervised_using_modified_label/Label_shuffled/output/BMMC_s2/meta_24.rds')

meta_BMMC_s_26 <- readRDS('./Benchmarking/Unsupervised/BMMC_s2/output/meta_26.rds')
unique(meta_BMMC_s_26$batch)

meta_BMMC_s_26$label_shuffled5 <- c(meta_c1$label_shuffled5,meta_c1$label_shuffled5,
                                meta_c2$label_shuffled5,
                                meta_c3$label_shuffled5,meta_c3$label_shuffled5,
                                meta_m1$label_shuffled5,meta_m1$label_shuffled5,
                                meta_m2$label_shuffled5,meta_m2$label_shuffled5,
                                meta_m3$label_shuffled5)

meta_BMMC_s_26$label_shuffled10 <- c(meta_c1$label_shuffled10,meta_c1$label_shuffled10,
                                meta_c2$label_shuffled10,
                                meta_c3$label_shuffled10,meta_c3$label_shuffled10,
                                meta_m1$label_shuffled10,meta_m1$label_shuffled10,
                                meta_m2$label_shuffled10,meta_m2$label_shuffled10,
                                meta_m3$label_shuffled10)

meta_BMMC_s_26$label_shuffled20 <- c(meta_c1$label_shuffled20,meta_c1$label_shuffled20,
                                meta_c2$label_shuffled20,
                                meta_c3$label_shuffled20,meta_c3$label_shuffled20,
                                meta_m1$label_shuffled20,meta_m1$label_shuffled20,
                                meta_m2$label_shuffled20,meta_m2$label_shuffled20,
                                meta_m3$label_shuffled20)
saveRDS(meta_BMMC_s_26,
        './Benchmarking/Supervised_using_modified_label/Label_shuffled/output/BMMC_s2/meta_26.rds')

meta_BMMC_s_35 <- readRDS('./Benchmarking/Unsupervised/BMMC_s2/output/meta_35.rds')
unique(meta_BMMC_s_35$batch)

meta_BMMC_s_35$label_shuffled5 <- c(meta_c1$label_shuffled5,meta_c1$label_shuffled5,
                                meta_c2$label_shuffled5,meta_c2$label_shuffled5,
                                meta_c3$label_shuffled5,
                                 
                                meta_m1$label_shuffled5,meta_m1$label_shuffled5,
                                meta_m2$label_shuffled5,
                                meta_m3$label_shuffled5,meta_m3$label_shuffled5)

meta_BMMC_s_35$label_shuffled10 <- c(meta_c1$label_shuffled10,meta_c1$label_shuffled10,
                                meta_c2$label_shuffled10,meta_c2$label_shuffled10,
                                meta_c3$label_shuffled10,
                                meta_m1$label_shuffled10,meta_m1$label_shuffled10,
                                meta_m2$label_shuffled10,
                                meta_m3$label_shuffled10,meta_m3$label_shuffled10)

meta_BMMC_s_35$label_shuffled20 <- c(meta_c1$label_shuffled20,meta_c1$label_shuffled20,
                                meta_c2$label_shuffled20,meta_c2$label_shuffled20,
                                meta_c3$label_shuffled20,
                                meta_m1$label_shuffled20,meta_m1$label_shuffled20,
                                meta_m2$label_shuffled20,
                                meta_m3$label_shuffled20,meta_m3$label_shuffled20)
saveRDS(meta_BMMC_s_35,
        './Benchmarking/Supervised_using_modified_label/Label_shuffled/output/BMMC_s2/meta_35.rds')

meta_BMMC_s_36 <- readRDS('./Benchmarking/Unsupervised/BMMC_s2/output/meta_36.rds')
unique(meta_BMMC_s_36$batch)

meta_BMMC_s_36$label_shuffled5 <- c(meta_c1$label_shuffled5,meta_c1$label_shuffled5,
                                meta_c2$label_shuffled5,meta_c2$label_shuffled5,
                                meta_c3$label_shuffled5,
                                 
                                meta_m1$label_shuffled5,meta_m1$label_shuffled5,
                                meta_m2$label_shuffled5,meta_m2$label_shuffled5,
                                meta_m3$label_shuffled5)

meta_BMMC_s_36$label_shuffled10 <- c(meta_c1$label_shuffled10,meta_c1$label_shuffled10,
                                meta_c2$label_shuffled10,meta_c2$label_shuffled10,
                                meta_c3$label_shuffled10,
                                meta_m1$label_shuffled10,meta_m1$label_shuffled10,
                                meta_m2$label_shuffled10,meta_m2$label_shuffled10,
                                meta_m3$label_shuffled10)

meta_BMMC_s_36$label_shuffled20 <- c(meta_c1$label_shuffled20,meta_c1$label_shuffled20,
                                meta_c2$label_shuffled20,meta_c2$label_shuffled20,
                                meta_c3$label_shuffled20,
                                meta_m1$label_shuffled20,meta_m1$label_shuffled20,
                                meta_m2$label_shuffled20,meta_m2$label_shuffled20,
                                meta_m3$label_shuffled20)
saveRDS(meta_BMMC_s_36,
        './Benchmarking/Supervised_using_modified_label/Label_shuffled/output/BMMC_s2/meta_36.rds')

# Retina

meta <- readRDS('./Benchmarking/Unsupervised/Retina/output/meta.rds')

head(meta)

meta$label_shuffled5 <- meta$celltype
meta$label_shuffled10 <- meta$celltype
meta$label_shuffled20 <- meta$celltype

sample_size <- ceiling(0.05 * nrow(meta))
set.seed(123)
idx <- sample(1:nrow(meta), sample_size,replace = F)
meta$label_shuffled5[idx] <- sample(meta$label_shuffled5[idx])

head(meta)

sample_size <- ceiling(0.1 * nrow(meta))
set.seed(123)
idx <- sample(1:nrow(meta), sample_size,replace = F)
meta$label_shuffled10[idx] <- sample(meta$label_shuffled10[idx])

sample_size <- ceiling(0.2 * nrow(meta))
set.seed(123)
idx <- sample(1:nrow(meta), sample_size,replace = F)
meta$label_shuffled20[idx] <- sample(meta$label_shuffled20[idx])

saveRDS(meta,'./Benchmarking/Supervised_using_modified_label/Label_shuffled/output/Retina/meta.rds')

meta_s1 <- readRDS('./Benchmarking/Unsupervised/Retina/output/meta_s1.rds')

all.equal(rownames(meta),rownames(meta_s1))

meta_s1$label_shuffled5 <- meta$label_shuffled5
meta_s1$label_shuffled10 <- meta$label_shuffled10
meta_s1$label_shuffled20 <- meta$label_shuffled20
saveRDS(meta_s1,
        './Benchmarking/Supervised_using_modified_label/Label_shuffled/output/Retina/meta_s1.rds')

meta_s2 <- readRDS('./Benchmarking/Unsupervised/Retina/output/meta_s2.rds')
all.equal(rownames(meta),rownames(meta_s2))
meta_s2$label_shuffled5 <- meta$label_shuffled5
meta_s2$label_shuffled10 <- meta$label_shuffled10
meta_s2$label_shuffled20 <- meta$label_shuffled20
saveRDS(meta_s2,
        './Benchmarking/Supervised_using_modified_label/Label_shuffled/output/Retina/meta_s2.rds')

meta_s3 <- readRDS('./Benchmarking/Unsupervised/Retina/output/meta_s3.rds')
all.equal(rownames(meta),rownames(meta_s3))
meta_s3$label_shuffled5 <- meta$label_shuffled5
meta_s3$label_shuffled10 <- meta$label_shuffled10
meta_s3$label_shuffled20 <- meta$label_shuffled20
saveRDS(meta_s3,
        './Benchmarking/Supervised_using_modified_label/Label_shuffled/output/Retina/meta_s3.rds')

meta_s4 <- readRDS('./Benchmarking/Unsupervised/Retina/output/meta_s4.rds')
all.equal(rownames(meta),rownames(meta_s4))
meta_s4$label_shuffled5 <- meta$label_shuffled5
meta_s4$label_shuffled10 <- meta$label_shuffled10
meta_s4$label_shuffled20 <- meta$label_shuffled20
saveRDS(meta_s4,
        './Benchmarking/Supervised_using_modified_label/Label_shuffled/output/Retina/meta_s4.rds')

# ABseq

meta <- readRDS('./Benchmarking/Unsupervised/Ab-seq/output/meta.rds')

head(meta)

meta$label_shuffled5 <- meta$celltype.l2
meta$label_shuffled10 <- meta$celltype.l2
meta$label_shuffled20 <- meta$celltype.l2

sample_size <- ceiling(0.05 * nrow(meta))
set.seed(123)
idx <- sample(1:nrow(meta), sample_size,replace = F)
meta$label_shuffled5[idx] <- sample(meta$label_shuffled5[idx])

sample_size <- ceiling(0.1 * nrow(meta))
set.seed(123)
idx <- sample(1:nrow(meta), sample_size,replace = F)
meta$label_shuffled10[idx] <- sample(meta$label_shuffled10[idx])

sample_size <- ceiling(0.2 * nrow(meta))
set.seed(123)
idx <- sample(1:nrow(meta), sample_size,replace = F)
meta$label_shuffled20[idx] <- sample(meta$label_shuffled20[idx])

saveRDS(meta,'./Benchmarking/Supervised_using_modified_label/Label_shuffled/output/Ab-seq/meta.rds')

meta_s1 <- readRDS('./Benchmarking/Unsupervised/Ab-seq/output/meta_s1.rds')
all.equal(rownames(meta),rownames(meta_s1))
meta_s1$label_shuffled5 <- meta$label_shuffled5
meta_s1$label_shuffled10 <- meta$label_shuffled10
meta_s1$label_shuffled20 <- meta$label_shuffled20
saveRDS(meta_s1,
        './Benchmarking/Supervised_using_modified_label/Label_shuffled/output/Ab-seq/meta_s1.rds')

meta_s2 <- readRDS('./Benchmarking/Unsupervised/Ab-seq/output/meta_s2.rds')
all.equal(rownames(meta),rownames(meta_s2))
meta_s2$label_shuffled5 <- meta$label_shuffled5
meta_s2$label_shuffled10 <- meta$label_shuffled10
meta_s2$label_shuffled20 <- meta$label_shuffled20
saveRDS(meta_s2,
        './Benchmarking/Supervised_using_modified_label/Label_shuffled/output/Ab-seq/meta_s2.rds')

meta_s3 <- readRDS('./Benchmarking/Unsupervised/Ab-seq/output/meta_s3.rds')
all.equal(rownames(meta),rownames(meta_s3))
meta_s3$label_shuffled5 <- meta$label_shuffled5
meta_s3$label_shuffled10 <- meta$label_shuffled10
meta_s3$label_shuffled20 <- meta$label_shuffled20
saveRDS(meta_s3,
        './Benchmarking/Supervised_using_modified_label/Label_shuffled/output/Ab-seq/meta_s3.rds')

meta_s4 <- readRDS('./Benchmarking/Unsupervised/Ab-seq/output/meta_s4.rds')
all.equal(rownames(meta),rownames(meta_s4))
meta_s4$label_shuffled5 <- meta$label_shuffled5
meta_s4$label_shuffled10 <- meta$label_shuffled10
meta_s4$label_shuffled20 <- meta$label_shuffled20
saveRDS(meta_s4,
        './Benchmarking/Supervised_using_modified_label/Label_shuffled/output/Ab-seq/meta_s4.rds')




meta <- readRDS('./Benchmarking/Unsupervised/TEA_s2/output/meta.rds')
head(meta)

meta$label_miss5 <- meta$celltype
meta$label_miss10 <- meta$celltype
meta$label_miss20 <- meta$celltype

sample_size <- ceiling(0.05 * nrow(meta))
set.seed(42)
idx <- sample(1:nrow(meta), sample_size,replace = F)
meta$label_miss5[idx] <- NA

sample_size <- ceiling(0.1 * nrow(meta))
set.seed(42)
idx <- sample(1:nrow(meta), sample_size,replace = F)
meta$label_miss10[idx] <- NA

sample_size <- ceiling(0.2 * nrow(meta))
set.seed(42)
idx <- sample(1:nrow(meta), sample_size,replace = F)
meta$label_miss20[idx] <- NA

saveRDS(meta,'./Benchmarking/Supervised_using_modified_label/Label_shuffled/output/TEA_s1/meta.rds')
saveRDS(meta,'./Benchmarking/Supervised_using_modified_label/Label_shuffled/output/TEA_s2/meta.rds')

meta <- readRDS('./Benchmarking/Reference-based_integration/output/Palette_Reference/meta.rds')
meta$label_miss5 <- meta$celltype.l2
meta$label_miss10 <- meta$celltype.l2
meta$label_miss20 <- meta$celltype.l2

sample_size <- ceiling(0.05 * nrow(meta))
set.seed(42)
idx <- sample(1:nrow(meta), sample_size,replace = F)
meta$label_miss5[idx] <- NA

sample_size <- ceiling(0.1 * nrow(meta))
set.seed(42)
idx <- sample(1:nrow(meta), sample_size,replace = F)
meta$label_miss10[idx] <- NA

sample_size <- ceiling(0.2 * nrow(meta))
set.seed(42)
idx <- sample(1:nrow(meta), sample_size,replace = F)
meta$label_miss20[idx] <- NA

meta_drop1 <- dplyr::filter(meta,batch != 'CITE1')
meta_drop2 <- dplyr::filter(meta,batch != 'CITE2')
meta_drop4 <- dplyr::filter(meta,batch != 'Multiome1')
meta_drop5 <- dplyr::filter(meta,batch != 'Multiome2')
meta_drop6 <- dplyr::filter(meta,batch != 'Multiome3')

saveRDS(meta_drop1,
        './Benchmarking/Supervised_using_modified_label/Label_missing/output/BMMC_s1/meta_drop1.rds')
saveRDS(meta_drop2,
        './Benchmarking/Supervised_using_modified_label/Label_missing/output/BMMC_s1/meta_drop2.rds')
saveRDS(meta_drop4,
        './Benchmarking/Supervised_using_modified_label/Label_missing/output/BMMC_s1/meta_drop4.rds')
saveRDS(meta_drop5,
        './Benchmarking/Supervised_using_modified_label/Label_missing/output/BMMC_s1/meta_drop5.rds')
saveRDS(meta_drop6,
        './Benchmarking/Supervised_using_modified_label/Label_missing/output/BMMC_s1/meta_drop6.rds')

meta_c1 <- dplyr::filter(meta,batch == 'CITE1')
meta_c2 <- dplyr::filter(meta,batch == 'CITE2')
meta_c3 <- dplyr::filter(meta,batch == 'CITE3')

meta_m1 <- dplyr::filter(meta,batch == 'Multiome1')
meta_m2 <- dplyr::filter(meta,batch == 'Multiome2')
meta_m3 <- dplyr::filter(meta,batch == 'Multiome3')

meta_BMMC_s <- readRDS('./Benchmarking/Unsupervised/BMMC_s2/output/meta_34.rds')

meta_BMMC_s$label_miss5 <- c(meta_c1$label_miss5,meta_c1$label_miss5,
                                meta_c2$label_miss5,meta_c2$label_miss5,
                                meta_c3$label_miss5,
                                meta_m1$label_miss5,
                                meta_m2$label_miss5,meta_m2$label_miss5,
                                meta_m3$label_miss5,meta_m3$label_miss5)

meta_BMMC_s$label_miss10 <- c(meta_c1$label_miss10,meta_c1$label_miss10,
                                meta_c2$label_miss10,meta_c2$label_miss10,
                                meta_c3$label_miss10,
                                meta_m1$label_miss10,
                                meta_m2$label_miss10,meta_m2$label_miss10,
                                meta_m3$label_miss10,meta_m3$label_miss10)

meta_BMMC_s$label_miss20 <- c(meta_c1$label_miss20,meta_c1$label_miss20,
                                meta_c2$label_miss20,meta_c2$label_miss20,
                                meta_c3$label_miss20,
                                meta_m1$label_miss20,
                                meta_m2$label_miss20,meta_m2$label_miss20,
                                meta_m3$label_miss20,meta_m3$label_miss20)

saveRDS(meta_BMMC_s,
        './Benchmarking/Supervised_using_modified_label/Label_missing/output/BMMC_s2/meta_34.rds')


meta_BMMC_s_24 <- readRDS('./Benchmarking/Unsupervised/BMMC_s2/output/meta_24.rds')

meta_BMMC_s_24$label_miss5 <- c(meta_c1$label_miss5,meta_c1$label_miss5,
                                meta_c2$label_miss5,
                                meta_c3$label_miss5,meta_c3$label_miss5,
                                meta_m1$label_miss5,
                                meta_m2$label_miss5,meta_m2$label_miss5,
                                meta_m3$label_miss5,meta_m3$label_miss5)

meta_BMMC_s_24$label_miss10 <- c(meta_c1$label_miss10,meta_c1$label_miss10,
                                meta_c2$label_miss10,
                                meta_c3$label_miss10,meta_c3$label_miss10,
                                meta_m1$label_miss10,
                                meta_m2$label_miss10,meta_m2$label_miss10,
                                meta_m3$label_miss10,meta_m3$label_miss10)

meta_BMMC_s_24$label_miss20 <- c(meta_c1$label_miss20,meta_c1$label_miss20,
                                meta_c2$label_miss20,
                                meta_c3$label_miss20,meta_c3$label_miss20,
                                meta_m1$label_miss20,
                                meta_m2$label_miss20,meta_m2$label_miss20,
                                meta_m3$label_miss20,meta_m3$label_miss20)
saveRDS(meta_BMMC_s_24,
        './Benchmarking/Supervised_using_modified_label/Label_missing/output/BMMC_s2/meta_24.rds')


meta_BMMC_s_26 <- readRDS('./Benchmarking/Unsupervised/BMMC_s2/output/meta_26.rds')
meta_BMMC_s_26$label_miss5 <- c(meta_c1$label_miss5,meta_c1$label_miss5,
                                meta_c2$label_miss5,
                                meta_c3$label_miss5,meta_c3$label_miss5,
                                meta_m1$label_miss5,meta_m1$label_miss5,
                                meta_m2$label_miss5,meta_m2$label_miss5,
                                meta_m3$label_miss5)

meta_BMMC_s_26$label_miss10 <- c(meta_c1$label_miss10,meta_c1$label_miss10,
                                meta_c2$label_miss10,
                                meta_c3$label_miss10,meta_c3$label_miss10,
                                meta_m1$label_miss10,meta_m1$label_miss10,
                                meta_m2$label_miss10,meta_m2$label_miss10,
                                meta_m3$label_miss10)

meta_BMMC_s_26$label_miss20 <- c(meta_c1$label_miss20,meta_c1$label_miss20,
                                meta_c2$label_miss20,
                                meta_c3$label_miss20,meta_c3$label_miss20,
                                meta_m1$label_miss20,meta_m1$label_miss20,
                                meta_m2$label_miss20,meta_m2$label_miss20,
                                meta_m3$label_miss20)
saveRDS(meta_BMMC_s_26,
        './Benchmarking/Supervised_using_modified_label/Label_missing/output/BMMC_s2/meta_26.rds')

meta_BMMC_s_35 <- readRDS('./Benchmarking/Unsupervised/BMMC_s2/output/meta_35.rds')
meta_BMMC_s_35$label_miss5 <- c(meta_c1$label_miss5,meta_c1$label_miss5,
                                meta_c2$label_miss5,meta_c2$label_miss5,
                                meta_c3$label_miss5,
                                 
                                meta_m1$label_miss5,meta_m1$label_miss5,
                                meta_m2$label_miss5,
                                meta_m3$label_miss5,meta_m3$label_miss5)

meta_BMMC_s_35$label_miss10 <- c(meta_c1$label_miss10,meta_c1$label_miss10,
                                meta_c2$label_miss10,meta_c2$label_miss10,
                                meta_c3$label_miss10,
                                meta_m1$label_miss10,meta_m1$label_miss10,
                                meta_m2$label_miss10,
                                meta_m3$label_miss10,meta_m3$label_miss10)

meta_BMMC_s_35$label_miss20 <- c(meta_c1$label_miss20,meta_c1$label_miss20,
                                meta_c2$label_miss20,meta_c2$label_miss20,
                                meta_c3$label_miss20,
                                meta_m1$label_miss20,meta_m1$label_miss20,
                                meta_m2$label_miss20,
                                meta_m3$label_miss20,meta_m3$label_miss20)
saveRDS(meta_BMMC_s_35,
        './Benchmarking/Supervised_using_modified_label/Label_missing/output/BMMC_s2/meta_35.rds')

meta_BMMC_s_36 <- readRDS('./Benchmarking/Unsupervised/BMMC_s2/output/meta_36.rds')
meta_BMMC_s_36$label_miss5 <- c(meta_c1$label_miss5,meta_c1$label_miss5,
                                meta_c2$label_miss5,meta_c2$label_miss5,
                                meta_c3$label_miss5,
                                 
                                meta_m1$label_miss5,meta_m1$label_miss5,
                                meta_m2$label_miss5,meta_m2$label_miss5,
                                meta_m3$label_miss5)

meta_BMMC_s_36$label_miss10 <- c(meta_c1$label_miss10,meta_c1$label_miss10,
                                meta_c2$label_miss10,meta_c2$label_miss10,
                                meta_c3$label_miss10,
                                meta_m1$label_miss10,meta_m1$label_miss10,
                                meta_m2$label_miss10,meta_m2$label_miss10,
                                meta_m3$label_miss10)

meta_BMMC_s_36$label_miss20 <- c(meta_c1$label_miss20,meta_c1$label_miss20,
                                meta_c2$label_miss20,meta_c2$label_miss20,
                                meta_c3$label_miss20,
                                meta_m1$label_miss20,meta_m1$label_miss20,
                                meta_m2$label_miss20,meta_m2$label_miss20,
                                meta_m3$label_miss20)
saveRDS(meta_BMMC_s_36,
        './Benchmarking/Supervised_using_modified_label/Label_missing/output/BMMC_s2/meta_36.rds')

meta <- readRDS('./Benchmarking/Unsupervised/Retina/output/meta.rds')

meta$label_miss5 <- meta$celltype
meta$label_miss10 <- meta$celltype
meta$label_miss20 <- meta$celltype

sample_size <- ceiling(0.05 * nrow(meta))
set.seed(42)
idx <- sample(1:nrow(meta), sample_size,replace = F)
meta$label_miss5[idx] <- NA

sample_size <- ceiling(0.1 * nrow(meta))
set.seed(42)
idx <- sample(1:nrow(meta), sample_size,replace = F)
meta$label_miss10[idx] <- NA

sample_size <- ceiling(0.2 * nrow(meta))
set.seed(42)
idx <- sample(1:nrow(meta), sample_size,replace = F)
meta$label_miss20[idx] <- NA
saveRDS(meta,'./Benchmarking/Supervised_using_modified_label/Label_missing/output/Retina/meta.rds')

meta_s1 <- readRDS('./Benchmarking/Unsupervised/Retina/output/meta_s1.rds')
all.equal(rownames(meta),rownames(meta_s1))
meta_s1$label_miss5 <- meta$label_miss5
meta_s1$label_miss10 <- meta$label_miss10
meta_s1$label_miss20 <- meta$label_miss20
saveRDS(meta_s1,
        './Benchmarking/Supervised_using_modified_label/Label_missing/output/Retina/meta_s1.rds')

meta_s2 <- readRDS('./Benchmarking/Unsupervised/Retina/output/meta_s2.rds')
all.equal(rownames(meta),rownames(meta_s2))
meta_s2$label_miss5 <- meta$label_miss5
meta_s2$label_miss10 <- meta$label_miss10
meta_s2$label_miss20 <- meta$label_miss20
saveRDS(meta_s2,
        './Benchmarking/Supervised_using_modified_label/Label_missing/output/Retina/meta_s2.rds')

meta_s3 <- readRDS('./Benchmarking/Unsupervised/Retina/output/meta_s3.rds')
all.equal(rownames(meta),rownames(meta_s3))
meta_s3$label_miss5 <- meta$label_miss5
meta_s3$label_miss10 <- meta$label_miss10
meta_s3$label_miss20 <- meta$label_miss20
saveRDS(meta_s3,
        './Benchmarking/Supervised_using_modified_label/Label_missing/output/Retina/meta_s3.rds')

meta_s4 <- readRDS('./Benchmarking/Unsupervised/Retina/output/meta_s4.rds')
all.equal(rownames(meta),rownames(meta_s4))
meta_s4$label_miss5 <- meta$label_miss5
meta_s4$label_miss10 <- meta$label_miss10
meta_s4$label_miss20 <- meta$label_miss20
saveRDS(meta_s4,
        './Benchmarking/Supervised_using_modified_label/Label_missing/output/Retina/meta_s4.rds')

meta <- readRDS('./Benchmarking/Unsupervised/Ab-seq/output/meta.rds')
meta$label_miss5 <- meta$celltype.l2
meta$label_miss10 <- meta$celltype.l2
meta$label_miss20 <- meta$celltype.l2

sample_size <- ceiling(0.05 * nrow(meta))
set.seed(42)
idx <- sample(1:nrow(meta), sample_size,replace = F)
meta$label_miss5[idx] <- NA

sample_size <- ceiling(0.1 * nrow(meta))
set.seed(42)
idx <- sample(1:nrow(meta), sample_size,replace = F)
meta$label_miss10[idx] <- NA

sample_size <- ceiling(0.2 * nrow(meta))
set.seed(42)
idx <- sample(1:nrow(meta), sample_size,replace = F)
meta$label_miss20[idx] <- NA

saveRDS(meta,'./Benchmarking/Supervised_using_modified_label/Label_missing/output/Ab-seq/meta.rds')


meta_s1 <- readRDS('./Benchmarking/Unsupervised/Ab-seq/output/meta_s1.rds')
all.equal(rownames(meta),rownames(meta_s1))
meta_s1$label_miss5 <- meta$label_miss5
meta_s1$label_miss10 <- meta$label_miss10
meta_s1$label_miss20 <- meta$label_miss20
saveRDS(meta_s1,
        './Benchmarking/Supervised_using_modified_label/Label_missing/output/Ab-seq/meta_s1.rds')


meta_s2 <- readRDS('./Benchmarking/Unsupervised/Ab-seq/output/meta_s2.rds')
all.equal(rownames(meta),rownames(meta_s2))
meta_s2$label_miss5 <- meta$label_miss5
meta_s2$label_miss10 <- meta$label_miss10
meta_s2$label_miss20 <- meta$label_miss20
saveRDS(meta_s2,
        './Benchmarking/Supervised_using_modified_label/Label_missing/output/Ab-seq/meta_s2.rds')

meta_s3 <- readRDS('./Benchmarking/Unsupervised/Ab-seq/output/meta_s3.rds')
all.equal(rownames(meta),rownames(meta_s3))
meta_s3$label_miss5 <- meta$label_miss5
meta_s3$label_miss10 <- meta$label_miss10
meta_s3$label_miss20 <- meta$label_miss20
saveRDS(meta_s3,
        './Benchmarking/Supervised_using_modified_label/Label_missing/output/Ab-seq/meta_s3.rds')

meta_s4 <- readRDS('./Benchmarking/Unsupervised/Ab-seq/output/meta_s4.rds')
all.equal(rownames(meta),rownames(meta_s4))
meta_s4$label_miss5 <- meta$label_miss5
meta_s4$label_miss10 <- meta$label_miss10
meta_s4$label_miss20 <- meta$label_miss20
saveRDS(meta_s4,
        './Benchmarking/Supervised_using_modified_label/Label_missing/output/Ab-seq/meta_s4.rds')

