library(future)
library(Signac)
library(Seurat)
library(GenomicRanges)
# library(parallel)
library(EnsDb.Hsapiens.v86)
library(data.table)
plan("multicore", workers = 10)
# options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM


peak_merge <- function(X,first = ":",second = "-",MAX = 10000,MIN = 20){
  
  splitS <- paste0('[',first,second,']')
  df <- list()
  for(i in 1:length(X)){
    tmp <- as.data.frame(do.call(rbind, lapply(X[[i]], function(x) {
      parts <- unlist(strsplit(x, splitS))
      return(parts)
    })))
    colnames(tmp) <- c("seqname","start","end")
    df[[i]] <- makeGRangesFromDataFrame(tmp)
  }
  
  g_uni <- reduce(Reduce(c,df))
  g_uni <- g_uni[!(g_uni@seqnames %in% c("chrX", "chrY"))]
  peakwidths <- width(g_uni)
  peaks <- g_uni[peakwidths  < MAX & peakwidths > MIN]
  return(peaks)
}

peak_merge2 <- function(X,MAX = 10000,MIN = 50){
  g_uni <- reduce(Reduce(c,X))
  peakwidths <- width(g_uni)
  peaks <- g_uni[peakwidths  < MAX & peakwidths > MIN]
  return(peaks)
}


new_count <- function(frag_path, cell_name = NULL, peak, min_cells = 5){
  
  frags <- Signac::CreateFragmentObject(
    path = frag_path,
    cells = cell_name
  )
  
  atac_count <- Signac::FeatureMatrix(
    fragments = frags,
    features = peak,
    cells = cell_name
  )
  
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotation) <- "UCSC"
  genome(annotation) <- "hg38"
  # create atac assay and add it to the object
  atac_assay <- CreateChromatinAssay(
    counts = atac_count,
    min.cells = min_cells,
    genome = 'hg38',
    fragments = frags,
    annotation = annotation
  )
  atac <- CreateSeuratObject(
    counts = atac_assay,
    assay = 'atac',
  )
  atac <- NucleosomeSignal(atac)
  atac <- TSSEnrichment(atac)
  
  return(atac)
}

new_count2 <- function(frag_path, cell_name = NULL, peak){
  
  frags <- Signac::CreateFragmentObject(
    path = frag_path,
    cells = cell_name
  )
  
  atac_count <- Signac::FeatureMatrix(
    fragments = frags,
    features = peak,
    cells = cell_name
  )
  
  return(atac_count)
}


Peak_from_frag <- function(path,cells = NULL,macs2.path){
  
  frags <- CreateFragmentObject(path, cells = cells)
  peaks <- CallPeaks(frags,macs2.path = macs2.path)
  # remove peaks on non-autosomes and in genomic blacklist regions
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  peaks <- peaks[!(peaks@seqnames %in% c("chrX", "chrY"))]
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
  
  return(peaks)
}

gen_rna <- function(data, min_cells = 3) {
  rna <- CreateSeuratObject(
    counts = data,
    min.cells = min_cells,
    assay = "rna"
  )
  rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")
  return(rna)
}

############################ plot
library(ggplot2)
library(cowplot)
library(patchwork)
require(RcppAnnoy)
library(RColorBrewer)
require(ggplot2)
library(ggrastr)
library(viridis)
require(grid)
require(ggpubr)
######################################################## Visualization using R ggplot2 package
ALL_COLORS = c("#ADFF2F", "#E19F20", "#DE8EE8", "#66C2A5", "#EF9563", "#FF7F00", "#86115A", "#A2CD5A", "#CDC0B0", "#FFFF33", "#7AC5CD", "#8B8B00", "#8B7D6B", "#586500", "#A6D854", "#CD950C",
               "#F781BF", "#FF664F", "#377EB8", "#8B1A1A", "#8A2BE2", "#CD1076", "#8DA0CB", "#EE6A50", "#E5C494", "#875300", "#E41A1C", "#FFD92F", "#8B008B", "#A65628", "#6CA6CD",
               "#BDB76B", "#4DAF4A", "#63BAAB", "#FFB6C1", "#8B6914", "#E78AC3", "#A52A2A", "#984EA3", "#FC8D62")#"#999999", "#DCDCDC", "#B3B3B3", "#474747", 
ALL_COLORS = c("#ADFF2F", "#E19F20", "#DE8EE8", "#66C2A5", "#EF9563", "#FF7F00", 
               "#86115A", "#A2CD5A", "#CDC0B0", "#FFFF33", "#7AC5CD", "#8B8B00", 
               "#8B7D6B", "#586500", "#CD950C", "#F781BF", "#FF664F", '#DDEFFF',
               "#377EB8", "#8B1A1A", "#8A2BE2", "#CD1076", "#8DA0CB", "#EE6A50", 
               "#E5C494", '#1CE6FF', '#FF34FF', '#FF4A46', '#008941', '#006FA6',
               '#A079BF', '#CC0744', '#C0B9B2', '#C2FF99', '#001E09', '#00489C',
               '#6F0062', '#0CBD66', '#EEC3FF', '#456D75', '#B77B68', '#7A87A1',
               '#788D66', '#885578', '#FAD09F', '#FF8A9A', '#D157A0', '#BEC459',
               '#456648', '#0086ED', '#886F4C', '#34362D', '#B4A8BD', '#00A6AA',
               '#452C2C', '#636375', '#A3C8C9', '#FF913F', '#938A81', '#575329',
               '#00FECF', '#B05B6F', '#8CD0FF', '#3B9700', '#A30059', '#FFDBE5', 
               '#7A4900', '#0000A6', '#63FFAC', '#B79762', '#004D43', '#8FB0FF', 
               '#997D87', '#5A0007', '#809693', '#6A3A4C', '#1B4400', '#4FC601',
               '#3B5DFF', '#4A3B53', '#FF2F80', '#61615A', '#BA0900', '#6B7900',
               '#00C2A0', '#FFAA92', '#FF90C9', '#B903AA', '#D16100', "#A6D854",
               '#000035', '#7B4F4B', '#A1C299', '#300018', '#0AA6D8', '#013349', 
               '#00846F', '#372101', '#FFB500', '#C2FFED', "#875300", "#E41A1C", 
               "#FFD92F", "#8B008B", "#A65628", "#6CA6CD", "#BDB76B", "#4DAF4A",
               "#63BAAB", "#FFB6C1", "#8B6914", "#E78AC3", "#A52A2A", "#984EA3", "#FC8D62")

color1 <- c("#FA8415","#2E9D32","#237AB6","#C63135","#9667B9","#5C5A9D","#e499a6","#a6bddb")

color2 <- c("#C9D7D9","#358587",
            "#EFC2B0","#BA2A21",
            "#CED7BE","#368B35")
col <- c("#ABDDDE","#FAD510","#C6CDF7","#F4B5BD","#0A9F9D","#FAEED1","#005295","#E6A0C4","#C52E19",
         "#9C9D9D")
col_3 <- c("#8b4513", "#6b8e23", "#483d8b", "#bc8f8f", "#008080", "#000080", "#daa520", "#8fbc8f", "#8b008b",
            "#b03060", "#ff0000", "#00ff00", "#00fa9a", "#8a2be2", "#dc143c", "#00ffff", "#00bfff", "#0000ff",
            "#adff2f", "#ff7f50", "#ff00ff", "#1e90ff", "#f0e68c", "#ffff54", "#add8e6", "#ff1493", "#ee82ee")

# col_13_<- c("#00ff00", "#ff0000", "#0000ff", "#ffff00", "#ff69b4", "#00ffff", "#ff00ff", "#008000", "#6495ed",  "#4b0082", "#eee8aa", "#2f4f4f", "#8b4513")
col_4 <- c("#8FC36D", "#f54646", "#4472c4", "#fff300", "#ff69b4", "#ff00ff", "#14e6e6", "#008000", "#82B4ed",  "#D4aaff", "#eee8aa", "#2f4f4f", "#ad6800")

col_46 <- c("#8FC36D", "#f54646", "#4472c4", "#fff300", "#ff69b4", "#ff00ff", "#14e6e6", "#008000", "#82B4ed",  "#D4aaff", "#eee8aa", "#2f4f4f", "#ad6800", "#0000ff")

# col_9  <- c("#ff4500", "#006400", "#0000ff", "#ffd700", "#ff1493", "#00ffff", "#4169e1", "#00ff00", "#bc8f8f")

# col_8 <- c("#00ff00", "#ff0000", "#0000ff", "#006400", "#c71585", "#00ffff", "#1e90ff", "#ffd700")
col_47 <- c("#8FC36D", "#f54646", "#4472c4", "#ff00ff", "#82B4ed", "#D4aaff", "#008000", "#fff300")

col_48  <- c("#8FC36D", "#82B4ed", "#D4aaff", "#0000ff", "#f54646")

col_49  <- c("#00ff00", "#ff0000", "#0000ff", "#87cefa")


# 渐变配色
col_5 <- c('#e7e1ef','#d4b9da','#c994c7','#df65b0','#e7298a','#ce1256','#980043','#67001f')
col_6 <- c('#e5f5f9','#ccece6','#99d8c9','#66c2a4','#41ae76','#238b45','#006d2c','#00441b')
col_8 <- c('#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b')
col_9 <- c('#e0f3db','#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe','#0868ac','#084081')
col_10 <- c('#fee8c8','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#b30000','#7f0000')
col_11 <- c('#ece7f2','#d0d1e6','#a6bddb','#74a9cf','#3690c0','#0570b0','#045a8d','#023858')
col_12 <- c('#ece2f0','#d0d1e6','#a6bddb','#67a9cf','#3690c0','#02818a','#016c59','#014636')
col_13 <- c('#fde0dd','#fcc5c0','#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177','#49006a')
col_14 <- c('#f7fcb9','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#006837','#004529')
col_15 <- c('#edf8b1','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494','#081d58')
col_17 <- c('#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506')
col_18 <- c('#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')

col_19 <- c('#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b')
col_20 <- c('#e5f5e0','#c7e9c0','#a1d99b','#74c476','#41ab5d','#238b45','#006d2c','#00441b')
col_21 <- c('#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#a63603','#7f2704')
col_22 <- c('#efedf5','#dadaeb','#bcbddc','#9e9ac8','#807dba','#6a51a3','#54278f','#3f007d')
col_23 <- c('#fee0d2','#fcbba1','#fc9272','#fb6a4a','#ef3b2c','#cb181d','#a50f15','#67000d')

# 过渡对比配色
col_24 <- c('#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#c7eae5','#80cdc1','#35978f','#01665e','#003c30')
col_25 <- c('#8e0152','#c51b7d','#de77ae','#f1b6da','#fde0ef','#e6f5d0','#b8e186','#7fbc41','#4d9221','#276419')
col_26 <- c('#40004b','#762a83','#9970ab','#c2a5cf','#e7d4e8','#d9f0d3','#a6dba0','#5aae61','#1b7837','#00441b')
col_27 <- c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')
col_28 <- c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')
col_29 <- c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695')
col_30 <- c('#a50026','#d73027','#f46d43','#fdae61','#fee08b','#d9ef8b','#a6d96a','#66bd63','#1a9850','#006837')
col_31 <- c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')

#对比色
# 8色
col_32 <- c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f','#bf5b17','#666666')
col_33 <- c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666')
col_34 <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00')
col_35 <- c('#fbb4ae','#b3cde3','#ccebc5','#decbe4','#fed9a6','#ffffcc','#e5d8bd','#fddaec')
col_36 <- c('#b3e2cd','#fdcdac','#cbd5e8','#f4cae4','#e6f5c9','#fff2ae','#f1e2cc','#cccccc')
col_37 <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf')
col_38 <- c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494','#b3b3b3')
col_39 <- c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5')


#9色
col_40 <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6')
col_41 <- c('#fbb4ae','#b3cde3','#ccebc5','#decbe4','#fed9a6','#ffffcc','#e5d8bd','#fddaec','#f2f2f2')
col_42 <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999')
col_43 <- c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9')

#12色
col_44 <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')
col_45 <- c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f')


visual_plot <- function(res_method, meta, method.use = "", dataset.use = "", output.dir = NULL, reduction = "umap", color_used = NULL, 
                        color_Cell_type = NULL, color_Batch = NULL, celltype_title = NULL) {
  library(ggplot2)
  library(cowplot)
  if (is.null(color_used) == T) {
    color_used = ALL_COLORS
    # c("#FFFF33", "#EE6A50", "#CD950C", "#A2CD5A", "#8B1A1A", "#6CA6CD", "#FFB6C1", "#E41A1C", "#CD1076", "#999999", "#8B008B", "#63BAAB", "#FF7F00",
    #              "#E19F20", "#BDB76B", "#984EA3", "#8A2BE2", "#586500", "#FF664F", "#DE8EE8", "#ADFF2F", "#8B8B00", "#875300", "#4DAF4A", "#F781BF", "#DCDCDC",
    #              "#A65628", "#8B7D6B", "#86115A", "#474747", "#EF9563", "#CDC0B0", "#A52A2A", "#8B6914", "#7AC5CD", "#377EB8", "#B3B3B3", "#E5C494", "#FFD92F",
    #              "#A6D854", "#E78AC3", "#8DA0CB", "#FC8D62", "#66C2A5")
  }
  if (is.null(color_Cell_type) == T) {
    color_Cell_type <- color_used
  }
  if (is.null(color_Batch) == T) {
    color_Batch <- rev(color_used)
  }
  if (reduction == "umap") {
    library(umap)
    out_umap <- uwot::umap(t(res_method), n_neighbors = 30, learning_rate = 0.5, init = "laplacian", 
                           metric = 'cosine', fast_sgd = FALSE, n_sgd_threads = 1,
                           min_dist = .1, n_threads = 4, ret_model = TRUE)
    umap_df<- as.data.frame(out_umap$embedding)
    umap_df$'batchlb' <- as.factor(meta$batchlb)
    if (is.numeric(meta$CellType) == TRUE) {
      umap_df$'CellType' <- meta$CellType
    } else {
      umap_df$'CellType' <- as.factor(meta$CellType)
    }
    rownames(umap_df) <- meta$cell
    colnames(umap_df) <- c('UMAP_1', 'UMAP_2', 'batchlb', 'CellType')
    
    if (is.numeric(umap_df$'CellType') == TRUE) {
      p01 <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, colour = batchlb)) + geom_point(size = 1, alpha = 0.6) + theme_bw() 
      p01 <- p01 + labs(x = NULL,y = NULL,title = method.use)
      p01 <- p01 + theme(legend.title = element_text(size=17), 
                         legend.key.size = unit(1.1, "cm"),
                         legend.key.width = unit(0.5,"cm"), 
                         legend.text = element_text(size=14), 
                         plot.title = element_text(color="black", size=20, hjust = 0.5),
                         panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")
                         , axis.text.x = element_blank(), axis.text.y = element_blank(),
                         axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
      library(RColorBrewer)
      p02 <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) + geom_point(size = 1, aes(colour = CellType), alpha = 0.6) + theme_bw() +
        scale_colour_gradientn(colours = topo.colors(10))
      p02 <- p02 + labs(x = NULL,y = NULL,title = celltype_title)###########################celltype = NULL
      p02 <- p02 + theme(legend.title = element_text(size=17), 
                         legend.key.size = unit(1.1, "cm"),
                         legend.key.width = unit(0.5,"cm"), 
                         legend.text = element_text(size=14), 
                         plot.title = element_text(color="black", size=20, hjust = 0.5),
                         panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")
                         , axis.text.x = element_blank(), axis.text.y = element_blank(),
                         axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
      p03 <- plot_grid(p01, p02, nrow = 2, ncol = 1)
      
      if (is.null(output.dir) == FALSE) {
        png(paste0(output.dir, method.use, "_", dataset.use, ".png"), width = 2*1000, height = 800, res = 2*72)
        print(plot_grid(p01, p02))
        dev.off()
        
        pdf(paste0(output.dir, method.use, "_", dataset.use, ".pdf"), width=15, height=7, paper='special')
        print(plot_grid(p01, p02))
        dev.off()
      }
      
      res_p <- list(p01, p02, p03)
    } else {
      p01 <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, colour = batchlb)) + geom_point(size = 1, alpha = 0.6) + theme_bw() +
        scale_color_manual(values = color_Batch, name = 'Batch')
      p01 <- p01 + labs(x = NULL,y = NULL,title = method.use)
      p01 <- p01 + theme(legend.title = element_text(size=17), 
                         legend.key.size = unit(1.1, "cm"),
                         legend.key.width = unit(0.5,"cm"), 
                         legend.text = element_text(size=14), 
                         plot.title = element_text(color="black", size=20, hjust = 0.5),
                         panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")
                         , axis.text.x = element_blank(), axis.text.y = element_blank(),
                         axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
      
      p02 <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, colour = CellType)) + geom_point(size = 1, alpha = 0.6) + theme_bw() +
        scale_color_manual(values = color_Cell_type, name = 'Cell type')
      p02 <- p02 + labs(x = NULL,y = NULL,title = celltype_title)
      p02 <- p02 + theme(legend.title = element_text(size=17), 
                         legend.key.size = unit(1.1, "cm"),
                         legend.key.width = unit(0.5,"cm"), 
                         legend.text = element_text(size=14), 
                         plot.title = element_text(color="black", size=20, hjust = 0.5),
                         panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")
                         , axis.text.x = element_blank(), axis.text.y = element_blank(),
                         axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
      p03 <- plot_grid(p01, p02, nrow = 2, ncol = 1)
      
      if (is.null(output.dir) == FALSE) {
        png(paste0(output.dir, method.use, "_", dataset.use, ".png"), width = 2*1000, height = 800, res = 2*72)
        print(plot_grid(p01, p02))
        dev.off()
        
        pdf(paste0(output.dir, method.use, "_", dataset.use, ".pdf"), width=15, height=7, paper='special')
        print(plot_grid(p01, p02))
        dev.off()
      }
      
      res_p <- list(p01, p02, p03)
    }
  } else if (reduction == "tsne") {
    library(Rtsne)
    out_tsne <- Rtsne(t(res_method))
    tsne_df<- as.data.frame(out_tsne$Y)
    tsne_df$'batchlb' <- as.factor(meta$batchlb)
    if (is.numeric(meta$CellType) == TRUE) {
      tsne_df$'CellType' <- meta$CellType
    } else {
      tsne_df$'CellType' <- as.factor(meta$CellType)
    }
    rownames(tsne_df) <- meta$cell
    colnames(tsne_df) <- c('TSNE_1', 'TSNE_2', 'batchlb', 'CellType')
    
    if (is.numeric(tsne_df$'CellType') == TRUE) {
      p01 <- ggplot(tsne_df, aes(x = TSNE_1, y = TSNE_2, colour = batchlb)) + geom_point(alpha = 0.6) +
        scale_color_manual(values = color_Batch, name = 'Batch')
      p01 <- p01 + labs(x = NULL,y = NULL,title = method.use)
      p01 <- p01 + theme(legend.title = element_text(size=17), 
                         legend.key.size = unit(1.1, "cm"),
                         legend.key.width = unit(0.5,"cm"), 
                         legend.text = element_text(size=14), 
                         plot.title = element_text(color="black", size=20, hjust = 0.5),
                         panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")
                         , axis.text.x = element_blank(), axis.text.y = element_blank(),
                         axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
      library(RColorBrewer)
      p02 <- ggplot(tsne_df, aes(x = TSNE_1, y = TSNE_2)) + geom_point(aes(colour = CellType), alpha = 0.6) + theme_bw()  +
        scale_color_manual(values = color_Cell_type, name = 'Cell type')
      p02 <- p02 + labs(x = NULL,y = NULL,title = celltype_title)
      p02 <- p02 + theme(legend.title = element_text(size=17), 
                         legend.key.size = unit(1.1, "cm"),
                         legend.key.width = unit(0.5,"cm"), 
                         legend.text = element_text(size=14), 
                         plot.title = element_text(color="black", size=20, hjust = 0.5),
                         panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")
                         , axis.text.x = element_blank(), axis.text.y = element_blank(),
                         axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
      p03 <- plot_grid(p01, p02, nrow = 2, ncol = 1)
      
      if (is.null(output.dir) == FALSE) {
        png(paste0(output.dir, method.use, "_", dataset.use, ".png"), width = 2*1000, height = 800, res = 2*72)
        print(plot_grid(p01, p02))
        dev.off()
        
        pdf(paste0(output.dir, method.use, "_", dataset.use, ".pdf"), width=15, height=7, paper='special')
        print(plot_grid(p01, p02))
        dev.off()
      }
      
      res_p <- list(p01, p02, p03)
    } else {
      p01 <- ggplot(tsne_df, aes(x = TSNE_1, y = TSNE_2, colour = batchlb)) + geom_point(alpha = 0.6) + theme_bw() +
        scale_color_manual(values = color_Batch, name = 'Batch')
      p01 <- p01 + labs(x = NULL,y = NULL,title = method.use)
      p01 <- p01 + theme(legend.title = element_text(size=17), 
                         legend.key.size = unit(1.1, "cm"),
                         legend.key.width = unit(0.5,"cm"), 
                         legend.text = element_text(size=14), 
                         plot.title = element_text(color="black", size=20, hjust = 0.5),
                         panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")
                         , axis.text.x = element_blank(), axis.text.y = element_blank(),
                         axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
      
      p02 <- ggplot(tsne_df, aes(x = TSNE_1, y = TSNE_2, colour = CellType)) + geom_point(alpha = 0.6) + theme_bw() +
        scale_color_manual(values = color_Cell_type, name = 'Cell type')
      p02 <- p02 + labs(x = NULL,y = NULL,title = celltype_title)
      p02 <- p02 + theme(legend.title = element_text(size=17), 
                         legend.key.size = unit(1.1, "cm"),
                         legend.key.width = unit(0.5,"cm"), 
                         legend.text = element_text(size=14), 
                         plot.title = element_text(color="black", size=20, hjust = 0.5),
                         panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")
                         , axis.text.x = element_blank(), axis.text.y = element_blank(),
                         axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
      p03 <- plot_grid(p01, p02, nrow = 2, ncol = 1)
      
      if (is.null(output.dir) == FALSE) {
        png(paste0(output.dir, method.use, "_", dataset.use, ".png"), width = 2*1000, height = 800, res = 2*72)
        print(plot_grid(p01, p02))
        dev.off()
        
        pdf(paste0(output.dir, method.use, "_", dataset.use, ".pdf"), width=15, height=7, paper='special')
        print(plot_grid(p01, p02))
        dev.off()
      }
      
      res_p <- list(p01, p02, p03)
    }
  }
  return(res_p)
}

dim.plot <- function (object, meta, colFactor = NULL, col.rev = F, title = NULL, ord = NULL,
                      genes = NULL, legend = TRUE, Colors = NULL, size = 0.5, Alpha = 0.8, 
                      plot.ncol = NULL, raw.count = FALSE, exp.range = NULL, exp.col = "firebrick2", 
                      dim1 = "UMAP1", dim2 = "UMAP2", dim.name = F,
                      label = FALSE, adjust.label = 0.25, label.font = 5) 
{
  X = Y = label0 = as.character()
  exp.col = as.character(exp.col)
  adjust0 = as.numeric(adjust.label)
  font0 = as.integer(label.font)
  dimReduce0 = object
  size0 = as.numeric(size)
  
  m0 = data.frame(X = dimReduce0[, 1], Y = dimReduce0[, 2])
  draw = NULL
  colFactor = as.character(colFactor)
  colData = meta
  m0$factor = as.factor(colData[, colFactor])
  if (!is.null(ord)) {
    m0$factor <- factor(m0$factor, levels=ord)
  }
  Cols = ALL_COLORS
  if (!is.null(Colors)) Cols = Colors
  # c("#FFFF33", "#EE6A50", "#CD950C", "#A2CD5A", "#8B1A1A", "#6CA6CD", "#FFB6C1", "#E41A1C", "#CD1076", "#999999", "#8B008B", "#63BAAB", "#FF7F00",
  #        "#E19F20", "#BDB76B", "#984EA3", "#8A2BE2", "#586500", "#FF664F", "#DE8EE8", "#ADFF2F", "#8B8B00", "#875300", "#4DAF4A", "#F781BF", "#DCDCDC",
  #        "#A65628", "#8B7D6B", "#86115A", "#474747", "#EF9563", "#CDC0B0", "#A52A2A", "#8B6914", "#7AC5CD", "#377EB8", "#B3B3B3", "#E5C494", "#FFD92F",
  #        "#A6D854", "#E78AC3", "#8DA0CB", "#FC8D62", "#66C2A5")
  if (col.rev) {
    Cols = rev(Cols)
  }
  Cols0 = Cols[1:length(levels(m0$factor))]
  if (!isTRUE(dim.name)) {
    dim1 = NULL; dim2 = NULL
  }
  
  if (length(levels(m0$factor)) <= 35) {
    if (isTRUE(label)) {
      m1 = data.frame(X = rep(0, length(levels(m0$factor))), 
                      Y = rep(0, length(levels(m0$factor))), label0 = levels(m0$factor))
      m1$X = aggregate(m0$X, list(m0$factor), mean)[,2]
      m1$Y = aggregate(m0$Y, list(m0$factor), mean)[,2]
      plot = ggplot(m0, aes(X, Y)) + 
        geom_point(aes(color = as.factor(factor)), alpha = Alpha, size = size0) + 
        geom_text(data = m1, aes(X, Y, label = label0), nudge_x = adjust0, nudge_y = adjust0, size = font0) + 
        scale_color_manual(values = Cols0) + 
        labs(x = dim1, y = dim2, color = colFactor, title = title) + 
        theme_bw(base_size = 15, base_line_size = 0) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              plot.title = element_text(color="black", size=16, hjust = 0.5),
              axis.text.x = element_blank(), axis.text.y = element_blank(),
              axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) + 
        guides(colour = guide_legend(override.aes = list(size = 4), ncol = ifelse(length(levels(m0$factor)) <= 15, 1, 2)))
    }
    else {
      plot = ggplot(m0, aes(X, Y)) + 
        geom_point(aes(color = as.factor(factor)), alpha = Alpha, size = size0) + 
        scale_color_manual(values = Cols0) + 
        labs(x = dim1, y = dim2, color = colFactor, title = title) + 
        theme_bw(base_size = 15, base_line_size = 0) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              plot.title = element_text(color="black", size=16, hjust = 0.5),
              axis.text.x = element_blank(), axis.text.y = element_blank(),
              axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) + 
        guides(colour = guide_legend(override.aes = list(size = 4), ncol = ifelse(length(levels(m0$factor)) <= 15, 1, 2)))
    }
  }
  else {
    if (isTRUE(label)) {
      m1 = data.frame(X = rep(0, length(levels(m0$factor))), 
                      Y = rep(0, length(levels(m0$factor))), label0 = levels(m0$factor))
      m1$X = aggregate(m0$X, list(m0$factor), mean)[,2]
      m1$Y = aggregate(m0$Y, list(m0$factor), mean)[,2]
      plot = ggplot(m0, aes(X, Y)) + 
        geom_point(aes(color = as.factor(factor)), 
                   alpha = Alpha, size = size0) + 
        geom_text(data = m1, aes(X, Y, label = label0), 
                  nudge_x = adjust0, nudge_y = adjust0, size = font0) + 
        labs(x = dim1, y = dim2, color = colFactor, title = title) + 
        theme_bw(base_size = 15, base_line_size = 0) + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              plot.title = element_text(color="black", size=16, hjust = 0.5)) + 
        guides(colour = guide_legend(override.aes = list(size = 4), ncol = ifelse(length(levels(m0$factor)) <= 15, 1, 2)))
    }
    else {
      plot = ggplot(m0, aes(X, Y)) + 
        geom_point(aes(color = as.factor(factor)),alpha = Alpha, size = size0) + 
        labs(x = dim1, y = dim2, color = colFactor, title = title) + 
        theme_bw(base_size = 15, base_line_size = 0) + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              plot.title = element_text(color="black", size=16, hjust = 0.5)) + 
        guides(colour = guide_legend(override.aes = list(size = 4), ncol = ifelse(length(levels(m0$factor)) <= 15, 1, 2)))
    }
  }
  if (!isTRUE(legend)) plot = plot + theme(legend.position="none")
  return(plot)
}

dim.plot0 <- function (object, meta, colFactor = NULL, col.rev = F, title = NULL,title_col = 'black',
                       ord = NULL, genes = NULL, legend = TRUE, Colors = NULL, size = 0.5, 
                       Alpha = 0.8, plot.ncol = NULL, raw.count = FALSE, exp.range = NULL, 
                       exp.col = "firebrick2", dim1 = "UMAP1", dim2 = "UMAP2", dim.name = F, 
                       label = FALSE, adjust.label = 0.25, label.font = 5, dpi = 750) 
{
  X = Y = label0 = as.character()
  exp.col = as.character(exp.col)
  adjust0 = as.numeric(adjust.label)
  font0 = as.integer(label.font)
  dimReduce0 = object
  size0 = as.numeric(size)
  m0 = data.frame(X = dimReduce0[, 1], Y = dimReduce0[, 2])
  draw = NULL
  colFactor = as.character(colFactor)
  colData = meta
  m0$factor = as.factor(colData[, colFactor])
  if (!is.null(ord)) {
    m0$factor <- factor(m0$factor, levels = ord)
  }

  Cols = ALL_COLORS
  if (!is.null(Colors)) Cols = Colors
  # c("#FFFF33", "#EE6A50", "#CD950C", "#A2CD5A", "#8B1A1A", 
  #        "#6CA6CD", "#FFB6C1", "#E41A1C", "#CD1076",  
  #        "#8B008B", "#63BAAB", "#FF7F00", "#E19F20", "#BDB76B", 
  #        "#984EA3", "#8A2BE2", "#586500", "#FF664F", "#DE8EE8", 
  #        "#ADFF2F", "#8B8B00", "#875300", "#4DAF4A", "#F781BF", 
  #        "#A65628", "#8B7D6B", "#86115A", 
  #        "#EF9563", "#CDC0B0", "#A52A2A", "#8B6914", "#7AC5CD", 
  #        "#377EB8", "#E5C494", "#FFD92F", "#A6D854", 
  #        "#E78AC3", "#8DA0CB", "#FC8D62", "#66C2A5")
  if (col.rev) {
    Cols = rev(Cols)
  }
  Cols0 = Cols[1:length(levels(m0$factor))]
  if (!isTRUE(dim.name)) {
    dim1 = NULL
    dim2 = NULL
  }
  if (length(levels(m0$factor)) <= 115 | !is.null(Colors)) {
    if (isTRUE(label)) {
      m1 = data.frame(X = rep(0, length(levels(m0$factor))), 
                      Y = rep(0, length(levels(m0$factor))), label0 = levels(m0$factor))
      m1$X = aggregate(m0$X, list(m0$factor), mean)[, 2]
      m1$Y = aggregate(m0$Y, list(m0$factor), mean)[, 2]
      plot = ggplot(m0, aes(X, Y)) + geom_point_rast(aes(color = as.factor(factor)), raster.dpi = dpi, stroke = 0.3, shape = 16,
                                                     alpha = Alpha, size = size0) + geom_text(data = m1, 
                                                                                              aes(X, Y, label = label0), nudge_x = adjust0, 
                                                                                              nudge_y = adjust0, size = font0) + scale_color_manual(values = Cols0) + 
        labs(x = dim1, y = dim2, color = colFactor, title = title) + 
        theme_void() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              plot.title = element_text(color = title_col, 
                                        size = 16, hjust = 0.5), axis.text.x = element_blank(), 
              axis.text.y = element_blank(), axis.ticks.x = element_blank(), 
              axis.ticks.y = element_blank()) + guides(colour = guide_legend(override.aes = list(size = 4), 
                                                                             ncol = ifelse(length(levels(m0$factor)) <= 15, 
                                                                                           1, 2)))
    }
    else {
      plot = ggplot(m0, aes(X, Y)) + geom_point_rast(aes(color = as.factor(factor)), raster.dpi = dpi, stroke = 0.3, shape = 16,
                                                     alpha = Alpha, size = size0) + scale_color_manual(values = Cols0) + 
        labs(x = dim1, y = dim2, color = colFactor, title = title) + 
        theme_void() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              plot.title = element_text(color = title_col, 
                                        size = 16, hjust = 0.5), axis.text.x = element_blank(), 
              axis.text.y = element_blank(), axis.ticks.x = element_blank(), 
              axis.ticks.y = element_blank()) + guides(colour = guide_legend(override.aes = list(size = 4), 
                                                                             ncol = ifelse(length(levels(m0$factor)) <= 15, 
                                                                                           1, 2)))
    }
  }
  else {
    if (isTRUE(label)) {
      m1 = data.frame(X = rep(0, length(levels(m0$factor))), 
                      Y = rep(0, length(levels(m0$factor))), label0 = levels(m0$factor))
      m1$X = aggregate(m0$X, list(m0$factor), mean)[, 2]
      m1$Y = aggregate(m0$Y, list(m0$factor), mean)[, 2]
      plot = ggplot(m0, aes(X, Y)) + geom_point_rast(aes(color = as.factor(factor)), raster.dpi = dpi, stroke = 0.3, shape = 16,
                                                     alpha = Alpha, size = size0) + geom_text(data = m1, 
                                                                                              aes(X, Y, label = label0), nudge_x = adjust0, 
                                                                                              nudge_y = adjust0, size = font0) + labs(x = dim1, 
                                                                                                                                      y = dim2, color = colFactor, title = title) + 
        theme_void() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              plot.title = element_text(color = title_col, 
                                        size = 16, hjust = 0.5)) + guides(colour = guide_legend(override.aes = list(size = 4), 
                                                                                                ncol = ifelse(length(levels(m0$factor)) <= 15, 
                                                                                                              1, 2)))
    }
    else {
      plot = ggplot(m0, aes(X, Y)) + geom_point_rast(aes(color = as.factor(factor)), raster.dpi = dpi, stroke = 0.3, shape = 16,
                                                     alpha = Alpha, size = size0) + labs(x = dim1, 
                                                                                         y = dim2, color = colFactor, title = title) + 
        theme_void() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              plot.title = element_text(color = title_col, 
                                        size = 16, hjust = 0.5)) + guides(colour = guide_legend(override.aes = list(size = 4), 
                                                                                                ncol = ifelse(length(levels(m0$factor)) <= 15, 
                                                                                                              1, 2)))
    }
  }
  if (!isTRUE(legend)) 
    plot = plot + theme(legend.position = "none")
  return(plot)
}

dim.plot1 <- function (object, meta, colFactor = NULL, col.rev = F, title = NULL,
                       ord = NULL, genes = NULL, legend = TRUE, Colors = NULL, size = 0.5, 
                       Alpha = 0.8, plot.ncol = NULL, raw.count = FALSE, exp.range = NULL, 
                       exp.col = "firebrick2", dim1 = "UMAP1", dim2 = "UMAP2", dim.name = F, 
                       label = FALSE, adjust.label = 0.25, label.font = 5) 
{
  X = Y = label0 = as.character()
  exp.col = as.character(exp.col)
  adjust0 = as.numeric(adjust.label)
  font0 = as.integer(label.font)
  dimReduce0 = object
  size0 = as.numeric(size)
  m0 = data.frame(X = dimReduce0[, 1], Y = dimReduce0[, 2])
  draw = NULL
  colFactor = as.character(colFactor)
  colData = meta
  m0$factor = as.factor(colData[, colFactor])
  if (!is.null(ord)) {
    m0$factor <- factor(m0$factor, levels = ord)
  }
  Cols = ALL_COLORS
  if (!is.null(Colors)) Cols = Colors
  # c("#FFFF33", "#EE6A50", "#CD950C", "#A2CD5A", "#8B1A1A", 
  #        "#6CA6CD", "#FFB6C1", "#E41A1C", "#CD1076",  
  #        "#8B008B", "#63BAAB", "#FF7F00", "#E19F20", "#BDB76B", 
  #        "#984EA3", "#8A2BE2", "#586500", "#FF664F", "#DE8EE8", 
  #        "#ADFF2F", "#8B8B00", "#875300", "#4DAF4A", "#F781BF", 
  #        "#A65628", "#8B7D6B", "#86115A", 
  #        "#EF9563", "#CDC0B0", "#A52A2A", "#8B6914", "#7AC5CD", 
  #        "#377EB8", "#E5C494", "#FFD92F", "#A6D854", 
  #        "#E78AC3", "#8DA0CB", "#FC8D62", "#66C2A5")
  if (col.rev) {
    Cols = rev(Cols)
  }
  Cols0 = Cols[1:length(levels(m0$factor))]
  if (!isTRUE(dim.name)) {
    dim1 = NULL
    dim2 = NULL
  }
  if (length(levels(m0$factor)) <= 115 | !is.null(Colors)) {
    if (isTRUE(label)) {
      m1 = data.frame(X = rep(0, length(levels(m0$factor))), 
                      Y = rep(0, length(levels(m0$factor))), label0 = levels(m0$factor))
      m1$X = aggregate(m0$X, list(m0$factor), mean)[, 2]
      m1$Y = aggregate(m0$Y, list(m0$factor), mean)[, 2]
      plot = ggplot(m0, aes(X, Y)) + geom_point(aes(color = as.factor(factor)),
                                                alpha = Alpha, size = size0) + geom_text(data = m1, 
                                                                                         aes(X, Y, label = label0), nudge_x = adjust0, 
                                                                                         nudge_y = adjust0, size = font0) + scale_color_manual(values = Cols0) + 
        labs(x = dim1, y = dim2, color = colFactor, title = title) + 
        theme_void() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              plot.title = element_text(color = "black", 
                                        size = 16, hjust = 0.5), axis.text.x = element_blank(), 
              axis.text.y = element_blank(), axis.ticks.x = element_blank(), 
              axis.ticks.y = element_blank()) + guides(colour = guide_legend(override.aes = list(size = 4), 
                                                                             ncol = ifelse(length(levels(m0$factor)) <= 15, 
                                                                                           1, 2)))
    }
    else {
      plot = ggplot(m0, aes(X, Y)) + geom_point(aes(color = as.factor(factor)),
                                                alpha = Alpha, size = size0) + scale_color_manual(values = Cols0) + 
        labs(x = dim1, y = dim2, color = colFactor, title = title) + 
        theme_void() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              plot.title = element_text(color = "black", 
                                        size = 16, hjust = 0.5), axis.text.x = element_blank(), 
              axis.text.y = element_blank(), axis.ticks.x = element_blank(), 
              axis.ticks.y = element_blank()) + guides(colour = guide_legend(override.aes = list(size = 4), 
                                                                             ncol = ifelse(length(levels(m0$factor)) <= 15, 
                                                                                           1, 2)))
    }
  }
  else {
    if (isTRUE(label)) {
      m1 = data.frame(X = rep(0, length(levels(m0$factor))), 
                      Y = rep(0, length(levels(m0$factor))), label0 = levels(m0$factor))
      m1$X = aggregate(m0$X, list(m0$factor), mean)[, 2]
      m1$Y = aggregate(m0$Y, list(m0$factor), mean)[, 2]
      plot = ggplot(m0, aes(X, Y)) + geom_point(aes(color = as.factor(factor)),
                                                alpha = Alpha, size = size0) + geom_text(data = m1, 
                                                                                         aes(X, Y, label = label0), nudge_x = adjust0, 
                                                                                         nudge_y = adjust0, size = font0) + labs(x = dim1, 
                                                                                                                                 y = dim2, color = colFactor, title = title) + 
        theme_void() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              plot.title = element_text(color = "black", 
                                        size = 16, hjust = 0.5)) + guides(colour = guide_legend(override.aes = list(size = 4), 
                                                                                                ncol = ifelse(length(levels(m0$factor)) <= 15, 
                                                                                                              1, 2)))
    }
    else {
      plot = ggplot(m0, aes(X, Y)) + geom_point(aes(color = as.factor(factor)),
                                                alpha = Alpha, size = size0) + labs(x = dim1, 
                                                                                    y = dim2, color = colFactor, title = title) + 
        theme_void() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              plot.title = element_text(color = "black", 
                                        size = 16, hjust = 0.5)) + guides(colour = guide_legend(override.aes = list(size = 4), 
                                                                                                ncol = ifelse(length(levels(m0$factor)) <= 15, 
                                                                                                              1, 2)))
    }
  }
  if (!isTRUE(legend)) 
    plot = plot + theme(legend.position = "none")
  return(plot)
}

dim.plot2 <- function (object, meta, colFactor = NULL, col.rev = F, title = NULL, labels = NULL, 
                       ord = NULL, genes = NULL, legend = TRUE, Colors = NULL, size = 0.5, 
                       Alpha = 0.8, plot.ncol = NULL, raw.count = FALSE, exp.range = NULL, 
                       exp.col = "firebrick2", dim1 = "UMAP1", dim2 = "UMAP2", dim.name = F, 
                       label = FALSE, adjust.label = 0.25, label.font = 5, dpi = 750) 
{
  X = Y = label0 = as.character()
  exp.col = as.character(exp.col)
  adjust0 = as.numeric(adjust.label)
  font0 = as.integer(label.font)
  dimReduce0 = object
  size0 = as.numeric(size)
  m0 = data.frame(X = dimReduce0[, 1], Y = dimReduce0[, 2])
  draw = NULL
  colFactor = as.character(colFactor)
  colData = meta
  m0$factor = as.factor(colData[, colFactor])
  if (!is.null(ord)) {
    m0$factor <- factor(m0$factor, levels = ord)
  }
  Cols = ALL_COLORS
  if (!is.null(Colors)) Cols = Colors
  # c("#FFFF33", "#EE6A50", "#CD950C", "#A2CD5A", "#8B1A1A", 
  #        "#6CA6CD", "#FFB6C1", "#E41A1C", "#CD1076",  
  #        "#8B008B", "#63BAAB", "#FF7F00", "#E19F20", "#BDB76B", 
  #        "#984EA3", "#8A2BE2", "#586500", "#FF664F", "#DE8EE8", 
  #        "#ADFF2F", "#8B8B00", "#875300", "#4DAF4A", "#F781BF", 
  #        "#A65628", "#8B7D6B", "#86115A", 
  #        "#EF9563", "#CDC0B0", "#A52A2A", "#8B6914", "#7AC5CD", 
  #        "#377EB8", "#E5C494", "#FFD92F", "#A6D854", 
  #        "#E78AC3", "#8DA0CB", "#FC8D62", "#66C2A5")
  if (col.rev) {
    Cols = rev(Cols)
  }
  Cols0 = Cols[1:length(levels(m0$factor))]
  if (!isTRUE(dim.name)) {
    dim1 = NULL
    dim2 = NULL
  }
  if (length(levels(m0$factor)) <= 115 | !is.null(Colors)) {
    if (isTRUE(label)) {
      m1 = data.frame(X = rep(0, length(levels(m0$factor))), 
                      Y = rep(0, length(levels(m0$factor))), label0 = levels(m0$factor))
      m1$X = aggregate(m0$X, list(m0$factor), mean)[, 2]
      m1$Y = aggregate(m0$Y, list(m0$factor), mean)[, 2]
      plot = ggplot(m0, aes(X, Y)) + geom_point_rast(aes(color = as.factor(factor)), raster.dpi = dpi, stroke = 0.3, shape = 16,
                                                     alpha = Alpha, size = size0) + geom_text(data = m1, 
                                                                                              aes(X, Y, label = label0), nudge_x = adjust0, 
                                                                                              nudge_y = adjust0, size = font0) + scale_color_manual(values = Cols0, labels = labels) + 
        labs(x = dim1, y = dim2, color = colFactor, title = title) + 
        theme_void() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              plot.title = element_text(color = "black", 
                                        size = 16, hjust = 0.5), axis.text.x = element_blank(), 
              axis.text.y = element_blank(), axis.ticks.x = element_blank(), 
              axis.ticks.y = element_blank()) + guides(colour = guide_legend(override.aes = list(size = 4), 
                                                                             ncol = ifelse(length(levels(m0$factor)) <= 15, 
                                                                                           1, 2)))
    }
    else {
      plot = ggplot(m0, aes(X, Y)) + geom_point_rast(aes(color = as.factor(factor)), raster.dpi = dpi, stroke = 0.3, shape = 16,
                                                     alpha = Alpha, size = size0) + scale_color_manual(values = Cols0, labels = labels) + 
        labs(x = dim1, y = dim2, color = colFactor, title = title) + 
        theme_void() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              plot.title = element_text(color = "black", 
                                        size = 16, hjust = 0.5), axis.text.x = element_blank(), 
              axis.text.y = element_blank(), axis.ticks.x = element_blank(), 
              axis.ticks.y = element_blank()) + guides(colour = guide_legend(override.aes = list(size = 4), 
                                                                             ncol = ifelse(length(levels(m0$factor)) <= 15, 
                                                                                           1, 2)))
    }
  }
  else {
    if (isTRUE(label)) {
      m1 = data.frame(X = rep(0, length(levels(m0$factor))), 
                      Y = rep(0, length(levels(m0$factor))), label0 = levels(m0$factor))
      m1$X = aggregate(m0$X, list(m0$factor), mean)[, 2]
      m1$Y = aggregate(m0$Y, list(m0$factor), mean)[, 2]
      plot = ggplot(m0, aes(X, Y)) + geom_point_rast(aes(color = as.factor(factor)), raster.dpi = dpi, stroke = 0.3, shape = 16,
                                                     alpha = Alpha, size = size0) + geom_text(data = m1, 
                                                                                              aes(X, Y, label = label0), nudge_x = adjust0, 
                                                                                              nudge_y = adjust0, size = font0) + labs(x = dim1, 
                                                                                                                                      y = dim2, color = colFactor, title = title) + 
        theme_void() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              plot.title = element_text(color = "black", 
                                        size = 16, hjust = 0.5)) + guides(colour = guide_legend(override.aes = list(size = 4), 
                                                                                                ncol = ifelse(length(levels(m0$factor)) <= 15, 
                                                                                                              1, 2)))
    }
    else {
      plot = ggplot(m0, aes(X, Y)) + geom_point_rast(aes(color = as.factor(factor)), raster.dpi = dpi, stroke = 0.3, shape = 16,
                                                     alpha = Alpha, size = size0) + labs(x = dim1, 
                                                                                         y = dim2, color = colFactor, title = title) + 
        theme_void() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              plot.title = element_text(color = "black", 
                                        size = 16, hjust = 0.5)) + guides(colour = guide_legend(override.aes = list(size = 4), 
                                                                                                ncol = ifelse(length(levels(m0$factor)) <= 15, 
                                                                                                              1, 2)))
    }
  }
  if (!isTRUE(legend)) 
    plot = plot + theme(legend.position = "none")
  return(plot)
}

plot_confidence <- function(object, confidence, 
                            title = NULL, 
                            is.legend = FALSE, legend.title = NULL,
                            size = 0.1, Alpha = 1, dpi = 450) {
  dimReduce0 = object
  m0 = data.frame(X = dimReduce0[, 1], Y = dimReduce0[, 2])
  m0$draw = as.numeric(confidence)
  fig.size(4,4)
  h1 = ggplot(m0, aes(X, Y)) + 
    geom_point_rast(aes(color = draw), raster.dpi = dpi, stroke = 0.3, shape = 16, alpha = Alpha, size = size) + 
    #scale_color_gradient(low = "gray90", high = "firebrick2") +
    scale_color_viridis(option="plasma", direction = -1, breaks=seq(0,1,1)) +
    labs(x = NULL, y = NULL, color = legend.title, title = title) + 
    theme_void() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          plot.title = element_text(color = "black", 
                                    size = 16, hjust = 0.5),
          legend.box.margin=margin(0,0,0,-5),
          #legend.key.size = unit(0.4,"cm"),#legend ??С
          legend.key.width = unit(0.25,"cm"),
          legend.key.height = unit(0.25,"cm"),
          legend.title=element_text(size=8,face="plain",color="black"),#legend title??С
          legend.text= element_text(size=7,face="plain",color="black"),#legend???ݴ?С
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(), axis.ticks.x = element_blank(), 
          axis.ticks.y = element_blank(), legend.position=c(0.8, 0.9), legend.direction="horizontal")
  h1
}

feature.plot <- function (object, expr.mat, genes = NULL, legend = TRUE, 
                          Colors = NULL, size = 0.5, Alpha = 0.8, plot.ncol = NULL, 
                          raw.count = FALSE, exp.range = NULL, exp.col = "firebrick2") 
{
  X = Y = as.character()
  exp.col = as.character(exp.col)
  size0 = as.numeric(size)
  
  dimReduce0 = object
  m0 = data.frame(X = dimReduce0[, 1], Y = dimReduce0[, 2])
  draw = NULL
  if (is.null(genes)) {
    stop("Please input gene symbol")
  }
  
  count0 = expr.mat
  genes = intersect(genes, rownames(count0))
  if (length(genes) == 1) {
    count0 = as.matrix(count0[genes, ])
    if (is.null(exp.range)) {
      m0$draw = as.numeric(count0)
    }
    else {
      m0$draw = as.numeric(as.matrix(count0))
      min0 = exp.range[1]
      max0 = exp.range[2]
      m0$draw[m0$draw < min0] = min0
      m0$draw[m0$draw > max0] = max0
    }
    if (isTRUE(legend)) {
      plot = ggplot(m0, aes(X, Y)) + 
        geom_point(aes(color = draw), alpha = Alpha, size = size0) + 
        scale_color_gradient(low = "gray90", high = exp.col) +
        labs(x = "UMAP1", y = "UMAP2", color = "Expression", title = genes) +
        theme_bw(base_size = 12, base_line_size = 0) + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              plot.title = element_text(color="black", size=16, hjust = 0.5))
    }
    else {
      plot = ggplot(m0, aes(X, Y)) + 
        geom_point(aes(color = draw), alpha = Alpha, size = size0) + 
        scale_color_gradient(low = "gray90", high = exp.col, guide = FALSE) + 
        labs(x = "UMAP1", y = "UMAP2", color = "Expression", title = genes) + 
        theme_bw(base_size = 12, base_line_size = 0) + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              plot.title = element_text(color="black", size=16, hjust = 0.5))
    }
    return(plot)
  }
  else {
    count0 = as.matrix(count0[genes, ])
    if (is.null(exp.range)) {
      count0 = count0
    }
    else {
      min0 = exp.range[1]
      max0 = exp.range[2]
      count0[count0 < min0] = min0
      count0[count0 > max0] = max0
    }
    if (isTRUE(legend)) {
      plot0 = lapply(genes, function(z) {
        ggplot(m0, aes(X, Y)) + 
          geom_point(aes(color = count0[z,]), alpha = Alpha, size = size0) + 
          scale_color_gradient(low = "gray90", high = exp.col) + 
          labs(x = "UMAP1", y = "UMAP2", color = "Expression", title = z) +
          theme_bw(base_size = 12, base_line_size = 0) + 
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                plot.title = element_text(color="black", size=16, hjust = 0.5))
      })
    }
    else {
      plot0 = lapply(genes, function(z) {
        ggplot(m0, aes(X, Y)) + 
          geom_point(aes(color = count0[z, ]), alpha = Alpha, size = size0) + 
          scale_color_gradient(low = "gray90", high = exp.col, guide = FALSE) +
          labs(x = "UMAP1", y = "UMAP2", color = "Expression", title = z) +
          theme_bw(base_size = 12, base_line_size = 0) + 
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                plot.title = element_text(color="black", size=16, hjust = 0.5))
      })
    }
    names(plot0) = genes
    if (is.null(plot.ncol)) {
      return(do.call(grid.arrange, c(plot0, ncol = ceiling(sqrt(length(genes))))))
    }
    else {
      return(do.call(grid.arrange, c(plot0, ncol = plot.ncol)))
    }
  }
}

Color0 <- function(n, col.rev = F) {
  cc <- c("#E41A1C", "#FF664F", "#00525B", "#377EB8", "#004B80", "#875300", "#E19F20", "#4DAF4A", "#005600", "#2A628F", 
          "#984EA3", "#DE8EE8", "#00C9BF", "#005B58", "#FF7F00", "#006B5E", "#63BAAB", "#FFFF33", "#586500", "#006B5F",
          "#A65628", "#EF9563", "#00C9AA", "#005B46", "#F781BF", "#86115A", "#007664", "#00C9B2", "#999999", "#474747", "#00754B", "#00C896")
  color_used <- c(
    "cyan4",      "skyblue3",   "darkolivegreen3",   "lightpink",   "darkmagenta",   "brown",   "blueviolet", "bisque4",  "deeppink3",       "darkkhaki",      
    "dodgerblue4",     "goldenrod4",            "gainsboro",       "firebrick4",      "cadetblue3",
    "greenyellow",     "gray6",           "coral2",                     "yellow4",         
    "darkgoldenrod3",  "navy",            "deepskyblue3","antiquewhite3"
  )
  col_f <- function(x) {
    cl <- rgb2hsv(col2rgb(x))
    hsv(h = cl[1,1], s = cl[2,1], v = cl[3,1])
  }
  cc <- c(cc, sapply(color_used, col_f))
  #scales::show_col(cc[order(cc, decreasing = T)])
  cc <- cc[order(cc, decreasing = T)][1:36]
  names(cc) <- NULL
  colors <- c()
  for(i in 1:6) {
    for(j in 1:6) {
      colors <- c(colors, cc[6*(j-1) + i])
    }
  }
  colors <- c(colors, rev(brewer.pal(8, "Set2")))
  Col = col2rgb(colors)
  dist.color = as.matrix(dist(t(Col)))
  diag(dist.color) = 1e10
  while(length(colors) > n) {
    minCol = apply(dist.color, 1, FUN = min)
    ids = which(minCol == min(minCol))[1]
    dist.color = dist.color[-ids, -ids]
    colors = colors[-ids]
  }
  if(col.rev) {
    colors <- rev(colors)
  }
  return(colors)
}

FeaturePlots<- function(obj, features, metadata_column, ...){
  all_cells<- colnames(obj)
  groups<- levels(obj@meta.data[, metadata_column])
  
  # the minimal and maximal of the value to make the legend scale the same. 
  minimal<- min(obj[['RNA']]@data[feature, ])
  maximal<- max(obj[['RNA']]@data[feature, ])
  ps<- list()
  for (group in groups) {
    subset_indx<- obj@meta.data[, metadata_column] == group
    subset_cells<- all_cells[subset_indx]
    p<- FeaturePlot(obj, features = feature, cells= subset_cells, ...) +
      scale_color_viridis_c(limits=c(minimal, maximal), direction = 1) +
      ggtitle(group) +
      theme(plot.title = element_text(size = 10, face = "bold"))
    ps[[group]]<- p
  }
  return(ps)
}
FeaturePlots<- function(obj, features, metadata_column, ...){
  all_cells<- colnames(obj)
  groups<- levels(obj@meta.data[, metadata_column])
  # the minimal and maximal of the value to make the legend scale the same. 
  ps<- list()
  for (feature in features) {
    minimal<- max(-3, min(obj[['RNA']]@data[feature, ]))
    maximal<- min(3, max(obj[['RNA']]@data[feature, ]))
    subset_indx<- obj@meta.data[, metadata_column] == feature
    subset_cells<- all_cells[subset_indx]
    p<- FeaturePlot(obj, features = feature, cells= subset_cells, ...) +
      scale_color_viridis_c(limits=c(minimal, maximal), direction = 1) +
      ggtitle(feature) +
      theme(plot.title = element_text(size = 10, face = "bold"))
    ps[[feature]]<- p
  }
  return(ps)
}

palette_colors <- function(alpha) {
  color = c('#ff6666', '#a30000', '#ff9d00', '#c3ff00', '#3a470e', '#34a321', 
            '#2198a3', '#001547', '#3400a3', '#9666ff', '#470e28', '#a3416d')
  color_rgb <- t(col2rgb(color)/255)
  grDevices::rgb(color_rgb[, 1], color_rgb[, 2], color_rgb[, 3], alpha = alpha)
}

# test functions
ff <- function(i, j, meta) {
  d1 <- as.matrix(S[[i]][[j]])
  d2 <- d1
  a1 <- which(d2>0, arr.ind = T)
  a2 <- a1
  a2 <- as.data.frame(a2)
  print(nrow(a2))
  a3 <- a2
  batches <- unique(meta$batchlb)
  a3[,1] <- as.character(meta$CellType[which(meta$batchlb == batches[i])][a3[,1]])
  a3[,2] <- as.character(meta$CellType[which(meta$batchlb == batches[j])][a3[,2]])
  print(length(which(a3[,1] == a3[, 2])))
  print(length(which(a3[,1] == a3[, 2]))/nrow(a2))
  print(summary(as.factor(a3[,1])))
  print(sum(summary(as.factor(a3[,1]))))
  print(summary(as.factor(a3[,2])))
  print(sum(summary(as.factor(a3[,2]))))
}

fff <- function(S, meta, r = 0) {
  if (!"dgCMatrix" %in% class(S)){
    S <- as(S, "dgCMatrix")
  }
  uS <- upper_tri(S@x, S@p, S@i, ncol(S))
  S <- sparseMatrix(i = as.vector(uS$i), j = as.vector(uS$j), x = as.vector(uS$x),
                    dims = c(ncol(S), ncol(S)), repr = "C")
  if (!"dgCMatrix" %in% class(S)){
    S <- as(S, "dgCMatrix")
  }
  if (r > 0) {
    S@x <- pmax(0, S@x - r)
  }
  a1 <- sparse_positive(S@x, S@p, S@i, ncol(S), length(which(S@x>0)))
  a1 <- as.data.frame(a1)
  print(nrow(a1))
  a1[,1] <- as.character(meta$CellType[a1[,1]])
  a1[,2] <- as.character(meta$CellType[a1[,2]])
  print(length(which(a1[,1] == a1[, 2])))
  print(length(which(a1[,1] == a1[, 2]))/nrow(a1))
  print(summary(as.factor(a1[,1])))
  print(sum(summary(as.factor(a1[,1]))))
  print(summary(as.factor(a1[,2])))
  print(sum(summary(as.factor(a1[,2]))))
}
ff0 <- function(S, meta) {
  
  for (i in 1:(length(S)-1)) {
    for (j in (i+1):length(S)) {
      d1 <- as.matrix(S[[i]][[j]])
      d2 <- d1
      a1 <- which(d2>0, arr.ind = T)
      
    }
  }
  
  a2 <- a1
  a2 <- as.data.frame(a2)
  print(nrow(a2))
  a3 <- a2
  batches <- unique(meta$batchlb)
  a3[,1] <- as.character(meta$CellType[which(meta$batchlb == batches[i])][a3[,1]])
  a3[,2] <- as.character(meta$CellType[which(meta$batchlb == batches[j])][a3[,2]])
  print(length(which(a3[,1] == a3[, 2])))
  print(length(which(a3[,1] == a3[, 2]))/nrow(a2))
  print(summary(as.factor(a3[,1])))
  print(sum(summary(as.factor(a3[,1]))))
  print(summary(as.factor(a3[,2])))
  print(sum(summary(as.factor(a3[,2]))))
  a1[,1] <- as.character(meta$CellType[a1[,1]])
  a1[,2] <- as.character(meta$CellType[a1[,2]])
  print(length(which(a1[,1] == a1[, 2])))
  print(length(which(a1[,1] == a1[, 2]))/nrow(a1))
  print(summary(as.factor(a1[,1])))
  print(sum(summary(as.factor(a1[,1]))))
  print(summary(as.factor(a1[,2])))
  print(sum(summary(as.factor(a1[,2]))))
}

library(gridExtra)
library(grid)
library(ggplot2)

grid_arrange_shared_legend <- function(plot_list, ncol = length(plot_list), nrow = 1, position = c("bottom", "right")) {
  
  plots <- plot_list
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}

grid_arrange_shared_multi_legends <- function(p_list1, p_list2, position = c("bottom", "right")) {
  plots1 <- p_list1
  plots2 <- p_list2
  
  position <- match.arg(position)
  
  g1 <- ggplotGrob(plots1[[1]] + theme(legend.position = position))$grobs
  legend1 <- g1[[which(sapply(g1, function(x) x$name) == "guide-box")]]
  lheight1 <- sum(legend1$height)
  lwidth1 <- sum(legend1$width)
  gl1 <- lapply(plots1, function(x) x + theme(legend.position="none"))
  gl1 <- c(gl1, ncol = length(gl1), nrow = 1)
  
  g2 <- ggplotGrob(plots2[[1]] + theme(legend.position = position))$grobs
  legend2 <- g2[[which(sapply(g2, function(x) x$name) == "guide-box")]]
  lheight2 <- sum(legend2$height)
  lwidth2 <- sum(legend2$width)
  gl2 <- lapply(plots2, function(x) x + theme(legend.position="none"))
  gl2 <- c(gl2, ncol = length(gl2), nrow = 1)
  
  g_blank <- ggplotGrob(p_list1[[1]] + theme(legend.position = NULL))$grobs
  g_blank <- g_blank[[1]]
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(
                       arrangeGrob(do.call(arrangeGrob, gl1),legend1,g_blank,ncol = 1,heights = unit.c(unit(1, "npc") - lheight2, lheight1, lheight2 - lheight1)),
                       arrangeGrob(do.call(arrangeGrob, gl2),legend2,ncol = 2,widths = unit.c(unit(1, "npc") - lheight2, lheight2)),
                       ncol = 2),
                     "right" = arrangeGrob(
                       arrangeGrob(do.call(arrangeGrob, gl1),legend1,g_blank,ncol = 3,widths = unit.c(unit(1, "npc") - lwidth2, lwidth1, lwidth2 - lwidth1)),
                       arrangeGrob(do.call(arrangeGrob, gl2),legend2,ncol = 2,widths = unit.c(unit(1, "npc") - lwidth2, lwidth2)),
                       nrow = 2)) %>% as.ggplot()
  combined          
}


.scatter.density.pc <- function(
    pcs, 
    group, # CellType
    x.name = "gcPC1",
    y.name = "gcPC2",
    color, 
    strokeSize, 
    pointSize, 
    strokeColor,
    alpha,
    title
){
  pList <- list()
  p <- ggplot(mapping = aes(
    x = pcs[, 1], 
    y = pcs[, 2], 
    fill = group)) +
    xlab(x.name) +
    ylab(y.name) +
    geom_point(
      aes(fill = group), 
      pch = 21, 
      color = strokeColor, 
      stroke = strokeSize, 
      size = pointSize,
      alpha = alpha) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    ggtitle(title) +
    theme(
      legend.position = "right",
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black", size = 0.5),
      plot.title=element_text(size = 15, hjust=0.5, margin = margin(0,0,0.1,0, unit = "cm")),
      legend.background = element_blank(),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.key = element_blank(),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14)) +
    guides(fill = guide_legend(override.aes = list(size = 4))) + 
    scale_fill_manual(name = group, values = color)
  le <- ggpubr::get_legend(p)
  p <- p + theme(legend.position = "none")
  xdens <- cowplot::axis_canvas(p, axis = "x")+
    geom_density(
      mapping = aes(
        x = pcs[,1], 
        fill = group),
      alpha = 0.7, 
      size = 0.2
    ) +
    theme(legend.position = "none") +
    scale_fill_manual(values = color)
  
  ydens <- cowplot::axis_canvas(
    p, 
    axis = "y", 
    coord_flip = TRUE) +
    geom_density(
      mapping = aes(
        x = pcs[,2],
        fill = group),
      alpha = 0.7,
      size = 0.2) +
    theme(legend.position = "none") +
    scale_fill_manual(name = group, values = color) +
    coord_flip()
  
  p1 <- insert_xaxis_grob(
    p,
    xdens,
    grid::unit(.2, "null"),
    position = "top"
  )
  p2 <- insert_yaxis_grob(
    p1,
    ydens,
    grid::unit(.2, "null"),
    position = "right"
  )
  pList[[1]] <- ggdraw(p2)
  pList[[2]] <- le
  return(pList)
}
