library(randomForest)
library(varImp)
library(caret)
library(Rcpp)

setwd('./MERFISH_ATAC')
#Random Forest Modelling
merf = readRDS("./MERFISH/merfish.rds")
rna = readRDS("./ATAC/rna.rds")
co_feat = intersect(rownames(merf),rownames(rna))
rm(rna)
merf = merf[co_feat,]
meta = readRDS("./All_anno.rds")
meta <- dplyr::filter(meta,modality == 'merfish')

all.equal(colnames(merf),rownames(meta))

merf_sub = as.data.frame(scale(t(as.matrix(merf))))
merf_sub$ann = as.factor(meta$Palette_anno) # for random forest prediction, using subclass
colnames(merf_sub) <- make.names(colnames(merf_sub))
model2 <- randomForest(ann ~ ., data = merf_sub, importance=TRUE)

i_scores <- caret::varImp(model2, conditional=TRUE) # calc the importance score of each marker
i_scores2 = MatrixGenerics::rowMaxs(as.matrix(i_scores))
df2 = data.frame(cbind(rownames(i_scores),i_scores2))
df2$i_scores2 = as.numeric(df2$i_scores2)

df_sorted <- df2[order(df2$i_scores2, decreasing = TRUE), ]
saveRDS(df_sorted,'./RF_all.rds')