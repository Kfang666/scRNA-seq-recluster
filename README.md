# scRecluster: Iterative Re-clustering for scRNA-seq Data
[![R-CMD-check](https://github.com/Kfang666/scRecluster/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Kfang666/scRecluster/actions)
[![License: GPL-3](https://img.shields.io/badge/License-GPL3-blue.svg)](https://opensource.org/licenses/GPL-3.0)
[![BioC status](http://www.bioconductor.org/shields/build/release/bioc/scRecluster.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/scRecluster)
> Hierarchical clustering refinement tool for single-cell RNA sequencing data
## 📦 Installation
### Development Version (GitHub)
```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("Kfang666/scRecluster")

#Key Features
Prevents over-clustering artifacts by automatically splitting oversized cell populations
Iterative hierarchical clustering: Splits clusters exceeding user-defined cell count (default: >2000 cells)
Seamless Seurat integration: Works directly with Seurat objects
Built-in visualization: One-command UMAP/tSNE cluster plotting
Marker gene export: Automatically generates CSV-formatted differential expression tables

#Quick Start
library(scRecluster)
library(Seurat)
# 1. Load data
data(pbmc3k, package = "Seurat")
pbmc <- pbmc3k
# 2. Run reclustering
pbmc <- iterative_clustering(
  object = pbmc,
  cell_threshold = 200,  # Split clusters >200 cells
  resolution = 0.6,
  pca_dims = 1:15
)
# 3. Visualize results
plot_clusters(pbmc, reduction = "umap")
# 4. Export marker genes
export_markers(pbmc, "markers.csv")


#Core Function Parameters
iterative_clustering(
  object,                 # Seurat object
  cluster_col = "seurat_clusters",  # Original cluster column
  cell_threshold = 2000,  # Maximum cells per cluster
  pca_dims = 1:30,        # PCA dimensions to use
  resolution = 1,         # Clustering resolution
  max_iterations = 1000,  # Maximum iterations
  parent_marker = "_",    # Hierarchy delimiter
  seed = 42               # Random seed
)
#Output Specifications
New metadata column: seurat_clusters_recluster (format: parentCluster_childCluster)
Visualization: plot_clusters() uses ggplot2 syntax for easy customization
Marker gene tables include: p_val, avg_log2FC, pct.1, pct.2, p_val_adj


#Dependencies
Imports:
    Seurat (>= 4.1.0),
    dplyr,
    ggplot2,
    utils

#Citation
@Manual{,
  title = {scRecluster: Hierarchical clustering refinement for single-cell RNA-seq},
  author = {Kai Fang},
  year = {2025},
  note = {R package version 1.0.0},
  url = {https://github.com/Kfang666/scRecluster},
}


