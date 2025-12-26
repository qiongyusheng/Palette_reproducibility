source("/home/server/sqy/code/midas/preprocess/utils.R")
# setwd("/root/workspace/code/midas/")
library(RColorBrewer)

parser <- ArgumentParser()
parser$add_argument("--task", type = "character", default = "TEA_s2_v3")
parser$add_argument("--method", type = "character", default = "midas_embed")
parser$add_argument("--experiment", type = "character", default = "e0")
parser$add_argument("--model", type = "character", default = "default")
parser$add_argument("--init_model", type = "character", default = "sp_latest")
o <- parser$parse_known_args()[[1]]

config <- parseTOML("/home/server/sqy/code/midas/configs/data.toml")[[gsub("_vd.*|_vt.*|_transfer$|_ref_.*$", "", o$task)]]
subset_names <- basename(config$raw_data_dirs)
subset_ids <- sapply(seq_along(subset_names) - 1, toString)
input_dirs <- pj("result", o$task, o$experiment, o$model, "predict", o$init_model, paste0("subset_", subset_ids))
pp_dir <- pj("data", "processed", o$task)
output_dir <- pj("result", "comparison", o$task, o$method, o$experiment, o$model, o$init_model)

mkdir(output_dir, remove_old = F)

label_paths <- pj(config$raw_data_dirs, "meta.rds")
K <- parseTOML("/home/server/sqy/code/midas/configs/model.toml")[["default"]]$dim_c
l <- 7.5  # figure size
L <- 10   # figure size
m <- 0.6  # legend margin

z_list <- list()
cell_name_list <- list()
label_list <- list()
is_label <- T
subset_name_list <- list()
S <- length(subset_names)

for (i in seq_along(subset_names)) {
  subset_name <- subset_names[i]
  z_dir    <- pj(input_dirs[i], "z", "joint")
  fnames <- dir(path = z_dir, pattern = ".csv$")
  fnames <- str_sort(fnames, decreasing = F)
  
  z_subset_list <- list()
  N <- length(fnames)
  for (n in seq_along(fnames)) {
    message(paste0("Loading Subset ", i, "/", S, ", File ", n, "/", N))
    z_subset_list[[n]] <- read.csv(file.path(z_dir, fnames[n]), header = F)
  }
  z_list[[subset_name]] <- bind_rows(z_subset_list)
  
  cell_name_list[[subset_name]] <- read.csv(pj(pp_dir, paste0("subset_", subset_ids[i]),
                                               "cell_names.csv"), header = T)[, 2]
  if (file.exists(label_paths[i])) {
    label_list[[subset_name]] <- readRDS(label_paths[i])[,2]
  } else {
    is_label <- F
  }
  
  subset_name_list[[subset_name]] <- rep(subset_name, length(cell_name_list[[subset_name]]))
}
rna <- t(data.matrix(bind_rows(z_list))) * 0  # pseudo rna counts
colnames(rna) <- do.call("c", unname(cell_name_list))
rownames(rna) <- paste0("rna-", seq_len(nrow(rna)))
obj <- CreateSeuratObject(counts = rna, assay = "rna")

z <- data.matrix(bind_rows(z_list))
c <- z[, 1:20]
colnames(c) <- paste0("c_", seq_len(ncol(c)))
rownames(c) <- colnames(obj)
obj[["c"]] <- CreateDimReducObject(embeddings = c, key = "c_", assay = "rna")

u <- z[, (K+1):(K+2)]
colnames(u) <- paste0("u_", seq_len(ncol(u)))
rownames(u) <- colnames(obj)
obj[["u"]] <- CreateDimReducObject(embeddings = u, key = "u_", assay = "rna")

obj@meta.data$batch <- factor(x = do.call("c", unname(subset_name_list)), levels = subset_names)
table(obj@meta.data$batch)[unique(obj@meta.data$batch)]
if (is_label) {
  obj@meta.data$l1 <- do.call("c", unname(label_list))
} else {
  obj@meta.data$l1 <- obj@meta.data$batch
}
write.csv(c,"./Benchmarking/Unsupervised/TEA_s2/output/midas_v3.csv",row.names=FALSE)
