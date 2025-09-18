source("./midas/preprocess/utils.R")

task <- "bench"
config <- parseTOML("./data.toml")[[gsub("_vd.*|_vt.*", "", task)]]
combs <- config$combs
comb_ratios <- config$comb_ratios
mods_ <- unique(unlist(combs))
mods <- vector()
for (mod in c("atac", "rna", "adt")) {
  if (mod %in% mods_) {
    mods <- c(mods, mod)
  }
}

input_dirs <- pj(config$raw_data_dirs)
task_dir <- pj("data", "processed", task)
mkdir(task_dir, remove_old = T)
output_feat_dir <- pj(task_dir, "feat")
mkdir(output_feat_dir, remove_old = T)

for (dataset_id in seq_along(input_dirs)) {
  subset_id <- toString(dataset_id - 1)
  output_dir <- pj(task_dir, paste0("subset_", subset_id))
  mkdir(pj(output_dir, "mat"), remove_old = T)
  mkdir(pj(output_dir, "mask"), remove_old = T)
}


merge_counts <- function(mod) {
  
  prt("Processing ", mod, " data ...")
  sc_list <- list()
  feat_list <- list()
  subset_id <- "0"
  for (dataset_id in seq_along(input_dirs)) {
    comb <- combs[[dataset_id]]
    comb_ratio <- comb_ratios[[dataset_id]]
    fp <- pj(input_dirs[dataset_id], paste0(config[[mod]][dataset_id]))
    if (file.exists(fp)) {
      prt("Loading ", fp, " ...")
      sc <- readRDS(fp)
      cell_num <- dim(sc)[2]
      end_ids <- round(cumsum(comb_ratio) / sum(comb_ratio) * cell_num)
      start_ids <- c(1, end_ids + 1)
      for (split_id in seq_along(comb)) {
        if (mod %in% comb[[split_id]]) {
          cell_names <- colnames(sc)[start_ids[split_id]:end_ids[split_id]]
          sc_list[[subset_id]] <- sc[,cell_names]
          feat_list[[subset_id]] <- rownames(sc_list[[subset_id]])
        }
        subset_id <- toString(strtoi(subset_id) + 1)
      }
    }
    else {
      subset_id <- toString(strtoi(subset_id) + length(comb))
    }
  }
  feat_union <- Reduce(union, feat_list)
  feat_dims[[mod]] <<- length(feat_union)
  
  write.csv(feat_union, file = pj(output_feat_dir, paste0("feat_names_", mod, ".csv")))
  
  # Get feature masks for each subset
  mask_list <- list()
  for (subset_id in names(sc_list)) {
    mask_list[[subset_id]] <- as.integer(feat_union %in% rownames(sc_list[[subset_id]]))
  }

  # Save subsets
  for (subset_id in names(sc_list)) {
    prt("Saving subset ", subset_id, " ...")
    output_dir <- pj(task_dir, paste0("subset_", subset_id))
    output_mat_dir <- pj(output_dir, "mat")
    
    mat <- t(data.matrix(sc_list[[subset_id]]))  # N * D
    if (grepl("_vd", task)) {
      mat <- as.matrix(downsampleMatrix(mat, prop = depth_factors[strtoi(subset_id) + 1], bycol=F))
    }
    if (grepl("_vt", task)) {
      mat <- mat[cell_mask_list[[subset_id]], ]
    }
    # Save count data
    write.csv(mat, file = pj(output_mat_dir, paste0(mod, ".csv")))
    # Save cell IDs
    write.csv(rownames(mat), file = pj(output_dir, "cell_names.csv"))
    
    output_mask_dir <- pj(output_dir, "mask")
    mask <- t(data.matrix(mask_list[[subset_id]]))  # 1 * D
    # Save the feature mask
    write.csv(mask, file = pj(output_mask_dir, paste0(mod, ".csv")))
  }
}

merge_frags <- function() {
  mod <- "atac"
  # Load different subsets
  prt("Processing ", mod, " data ...")
  sc_list <- list()
  feat_list <- list()
  
  subset_id <- "0"
  for (dataset_id in seq_along(input_dirs)) {
    comb <- combs[[dataset_id]]
    comb_ratio <- comb_ratios[[dataset_id]]
    fp <- pj(input_dirs[dataset_id], paste0(config[[mod]][dataset_id]))
    if (file.exists(fp)) {
      prt("Loading ", fp, " ...")
      sc <- readRDS(fp)
      cell_num <- dim(sc)[2]
      end_ids <- round(cumsum(comb_ratio) / sum(comb_ratio) * cell_num)
      start_ids <- c(1, end_ids + 1)
      for (split_id in seq_along(comb)) {
        if (mod %in% comb[[split_id]]) {
          cell_names <- colnames(sc)[start_ids[split_id]:end_ids[split_id]]
          sc_list[[subset_id]] <- sc[,cell_names]# subset(sc, cells = cell_names)
          feat_list[[subset_id]] <- rownames(sc_list[[subset_id]])
        }
        subset_id <- toString(strtoi(subset_id) + 1)
      }
    }
    else {
      subset_id <- toString(strtoi(subset_id) + length(comb))
    }
  }
  feat_merged <- Reduce(union, feat_list) 

  gr_sorted <- sort(StringToGRanges(feat_merged))
  feat_dims[[mod]] <<- width(gr_sorted@seqnames)
  feat_sorted <- GRangesToString(gr_sorted)
  write.csv(feat_sorted, file = pj(output_feat_dir, paste0("feat_names_", mod, ".csv")))
  
  for (subset_id in names(sc_list)) {
    prt("Saving subset ", subset_id, " ...")
    output_dir <- pj(task_dir, paste0("subset_", subset_id))
    output_mat_dir <- pj(output_dir, "mat")
    
    mat <- t(data.matrix(sc_list[[subset_id]])[feat_merged, ])  # N * D
    if (grepl("_vd", task)) {
      mat <- as.matrix(downsampleMatrix(mat, prop = depth_factors[strtoi(subset_id) + 1], bycol=F))
    }
    if (grepl("_vt", task)) {
      mat <- mat[cell_mask_list[[subset_id]], ]
    }
    # Save count data
    write.csv(mat, file = pj(output_mat_dir, paste0(mod, ".csv")))
    # Save cell IDs
    write.csv(rownames(mat), file = pj(output_dir, "cell_names.csv"))
  }
}


feat_dims <- list()
for (mod in mods) {
  if (mod == "atac") {
    merge_frags()
  } else {
    merge_counts(mod)
  }
}

# Save feature dimensionalities
prt("feat_dims: ", feat_dims)
write.csv(feat_dims, file = pj(output_feat_dir, "feat_dims.csv"))
