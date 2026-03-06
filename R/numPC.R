#' Best Number of PCs
#'
#' This function runs PCA and determines the best number of principal components to use for finding neighbors, tSNE or UMAP analysis. It returns a number and plots the variance vs cumulative variance of each PC.
#' @param obj_Seurat Seurat object.
#' @param pca_name Name of PCA reduction to determine best number of PCs.
#' @keywords pca
#' @export
#'

numPC <- function(obj_Seurat, pca_name){
  # Determine percent of variation associated with each PC
  pct <- obj_Seurat@reductions[[pca_name]]@stdev^2 / sum(obj_Seurat@reductions[[pca_name]]@stdev^2) * 100

  # Calculate cumulative percents for each PC
  cum <- cumsum(pct)

  # Determine which PC exhibits cumulative percent greater than 80%
  numPC <- which(cum > 80 )[1]

  summary.pca <- data.frame(prop.var=pct, cum.prop=cum, pc=c(1:length(pct)))

  # Plot variance vs cumulative variance of each PC (first)
  #pdf(paste0(pca_name, "_cumulative_variance.pdf"))
  p <- ggplot(data=summary.pca, aes(x=pc)) + geom_col(aes(y=prop.var), fill = "blue") + geom_line(aes(y=cum.prop), color = "red") +
    geom_vline(aes(xintercept=numPC), linetype="dashed") + theme_bw()
  ggsave(paste0(pca_name, "_cumulative_variance.pdf"), plot = p)
  #invisible(dev.off())

  numPC
}
