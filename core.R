#' Iterative clustering refinement
#' 
#' @param object Seurat对象
#' @param cluster_col 原始聚类列名（默认："seurat_clusters"）
#' @param cell_threshold 单群最大细胞数（默认：2000）
#' @param pca_dims PCA使用维度（默认：1:30）
#' @param resolution 聚类分辨率（默认：1）
#' @param max_iterations 最大迭代次数（默认：1000）
#' @param parent_marker 层级标签分隔符（默认："_"）
#' @param seed 随机种子（默认：42）
#' @return 修改后的Seurat对象，包含新的seurat_clusters_recluster列
#' @export
#' @examples
#' \dontrun{
#' pbmc <- iterative_clustering(pbmc, cell_threshold = 1500)
#' }
iterative_clustering <- function(object,
                                 cluster_col = "seurat_clusters",
                                 cell_threshold = 2000,
                                 pca_dims = 1:30,
                                 resolution = 1,
                                 max_iterations = 1000,
                                 parent_marker = "_",
                                 seed = 42) {
  
  # 参数验证
  validate_params(cell_threshold, pca_dims, resolution, max_iterations)
  
  # 初始化结果列
  object$seurat_clusters_recluster <- as.character(object[[cluster_col, drop = TRUE]])
  
  # 主迭代循环
  for (iter in 1:max_iterations) {
    cluster_counts <- table(object$seurat_clusters_recluster)
    over_clusters <- cluster_counts[cluster_counts > cell_threshold]
    
    if (length(over_clusters) == 0) break
    
    object <- process_single_cluster(
      object = object,
      cluster_name = names(which.max(over_clusters)),
      pca_dims = pca_dims,
      resolution = resolution,
      parent_marker = parent_marker,
      seed = seed
    )
    
    message(sprintf("[Iteration %d] Current max cluster size: %d", 
                    iter, max(cluster_counts)))
  }
  
  return(object)
}

# 内部函数：参数验证
validate_params <- function(cell_threshold, pca_dims, resolution, max_iterations) {
  if (!is.numeric(cell_threshold)) stop("cell_threshold must be numeric")
  if (any(pca_dims < 1)) stop("PCA dimensions must be >=1")
  if (resolution <= 0) stop("Resolution must be positive")
  if (max_iterations < 1) stop("Max iterations must be >=1")
}

# 内部函数：处理单个聚类
process_single_cluster <- function(object, cluster_name, pca_dims, resolution, parent_marker, seed) {
  set.seed(seed)
  cells <- colnames(object)[object$seurat_clusters_recluster == cluster_name]
  sub_obj <- subset(object, cells = cells)
  
  Seurat::VariableFeatures(sub_obj) <- Seurat::VariableFeatures(object)
  
  sub_obj <- sub_obj %>%
    Seurat::DietSeurat(assays = "RNA") %>%
    Seurat::ScaleData(verbose = FALSE) %>%
    Seurat::RunPCA(features = Seurat::VariableFeatures(object), 
                   npcs = max(pca_dims),
                   verbose = FALSE) %>%
    Seurat::FindNeighbors(dims = pca_dims, verbose = FALSE) %>%
    Seurat::FindClusters(resolution = resolution, 
                         algorithm = 2,
                         verbose = FALSE)
  
  new_labels <- paste(cluster_name, sub_obj$seurat_clusters, sep = parent_marker)
  object$seurat_clusters_recluster[cells] <- new_labels
  
  return(object)
}