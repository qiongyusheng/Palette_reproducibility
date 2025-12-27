# Palette Reproducibility
Code and data to reproduce the results in Palette manuscript. Note in the manuscript we used a development version of Palette with slightly different grammar and can also be found there.

We also provide an R package called `PaletteSC` and detailed tutorials (available on [GitHub](https://github.com/qiongyusheng/Palette)).

## Code structure

  - `Palette`, `midas`, and `scVAEIT`
    contain the implementations of the three methods evaluated in our manuscript.
    
  - `metrics`
    contains the implementation of evaluation metrics used in the study.

  - `Benchmarking`
    contains code for comprehensive benchmarking experiments, including:
     1. Unsupervised mosaic integration
         - [Balanced cell type composition across batches](https://github.com/qiongyusheng/Palette_reproducibility/tree/main/Benchmarking/Unsupervised)
         - [Imbalanced cell type composition across batches](https://github.com/qiongyusheng/Palette_reproducibility/tree/main/Benchmarking/CellType_missing)
      2. Supervised mosaic integration
         - [Ground-truth labels](https://github.com/qiongyusheng/Palette_reproducibility/tree/main/Benchmarking/Supervised) 
         - [Low-quality labels](https://github.com/qiongyusheng/Palette_reproducibility/tree/main/Benchmarking/Supervised_using_modified_label/Label_shuffled)
         - [Partially missing labels](https://github.com/qiongyusheng/Palette_reproducibility/tree/main/Benchmarking/Supervised_using_modified_label/Label_miss) 
         - [Unimodal-derived labels](https://github.com/qiongyusheng/Palette_reproducibility/tree/main/Benchmarking/Supervised_using_modified_label/BMMC_s1_unimodal_label)
      3. [Reference-based integration](https://github.com/qiongyusheng/Palette_reproducibility/tree/main/Benchmarking/Reference-based_integration)

    
  - `Cross_condition_PBMC`
    contains code for the analysis of cross-condition human PBMC mosaic datasets.

  - `Cross_species`
    contains code for cross-species integration experiments, including  MOp multimodal and WAT scRNA-seq datasets.

  - `diagonal_integration`
    contains code for diagonal integration experiments using multiple datasets.

  - `Human_tonsil_10xVisium`
    contains code for the analysis of the human tonsil 10x Visium spatial dataset.
    
  - `scRNA_Batch_integration`
    contains code for unimodal scRNA-seq batch integration experiments using multiple datastes, including both unsupervised and supervised (ground-truth labels, low-quality labels, and partially missing labels) mode.

  - `modality_inference`
    contains code for imputing missing modalities in multiple multimodal single-cell datasets.

---

## Data Availability

 Processed datasets used for the analyses have been deposited to [Zenodo.](https://zenodo.org/records/18045027)
