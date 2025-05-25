#' Export marker genes
#' 
#' @param object Seurat对象
#' @param output_file 输出文件名（默认："markers.csv"）
#' @param min.pct 最小表达比例（默认：0.25）
#' @param logfc.threshold 最小logFC阈值（默认：0.585）
#' @param ... 传递给FindAllMarkers的额外参数
#' @export
export_markers <- function(object, 
                           output_file = "markers.csv",
                           min.pct = 0.25,
                           logfc.threshold = 0.585,
                           ...) {
  Seurat::Idents(object) <- "seurat_clusters_recluster"
  markers <- Seurat::FindAllMarkers(
    object,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold,
    ...
  )
  utils::write.csv(markers, file = output_file)
  return(invisible(markers))
}

#' Visualization of clustering results
#' 
#' @param object Seurat对象
#' @param ... 传递给DimPlot的额外参数
#' @export
plot_clusters <- function(object, ...) {
  Seurat::DimPlot(object, 
                  group.by = "seurat_clusters_recluster",
                  label = TRUE,
                  ...) + 
    ggplot2::ggtitle("Refined Clusters") +
    Seurat::NoLegend()
}