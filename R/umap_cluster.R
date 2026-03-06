#' Clustering and UMAP
#'
#' This function uses the output of seurat_PCA. It returns a Seurat object and can plot the clustering in UMAP.
#' @param pca_Seurat Seurat object with PCA reduction data.
#' @param assay Assay to perform analysis on.
#' @param pca_name Name of PCA data reduction to use in analysis.
#' @param umap_name Name of UMAP data reduction to be stored in the Seurat object.
#' @param resolution Resolution setting. Higher number is used for larger data sets (more clusters). Default 0.2.
#' @param plots Logical. Return UMAP clustering saved as PDF. Default to TRUE.
#' @keywords umap
#' @export
#'

umap_cluster <- function(pca_Seurat, assay, pca_name, umap_name, resolution=0.2, plots=TRUE){

  graph.name <- paste0(assay, "_snn")

  numPC <- numPC(pca_Seurat, pca_name)

  pca_Seurat <- FindNeighbors(pca_Seurat, assay = assay, reduction = pca_name, dims = 1:numPC)
  pca_Seurat <- FindClusters(pca_Seurat, resolution = resolution, graph.name = graph.name)
  pca_Seurat <- RunUMAP(pca_Seurat, dims = 1:numPC, reduction = pca_name, assay = assay, reduction.name = umap_name, verbose = FALSE)

  if(plots){
    #pdf(paste0(umap_name, "_umap.pdf"))
    DefaultAssay(pca_Seurat) <- assay
    p <- DimPlot(pca_Seurat, reduction = umap_name, label = TRUE, label.size = 6)
    ggsave(paste0(umap_name, "_umap.pdf"), plot = p)
    #invisible(dev.off())
  }

  pca_Seurat
}
