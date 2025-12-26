meta <- readRDS('./meta.rds')

meta$label_shuffled5 <- meta$CellType
meta$label_shuffled10 <- meta$CellType
meta$label_shuffled20 <- meta$CellType

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

meta$label_miss5 <- meta$CellType
meta$label_miss10 <- meta$CellType
meta$label_miss20 <- meta$CellType

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

meta$label_miss5s <- meta$label_miss5
meta$label_miss10s <- meta$label_miss10
meta$label_miss20s <- meta$label_miss20
meta$label_miss5s[which(is.na(meta$label_miss5s))] <- "Unknown"
meta$label_miss10s[which(is.na(meta$label_miss10s))] <- "Unknown"
meta$label_miss20s[which(is.na(meta$label_miss20s))] <- "Unknown"

saveRDS(meta,"./meta.rds")
