#' Scatter Plots of ADTs
#'
#' This function creates scatter plots for each combination of ADTs. It returns a directory containing all the plots.
#' @param obj_Seurat Seurat object with ADT assay.
#' @param no_ADTs Vector of ADTs to exclude from the graphs. Default will return graphs for all ADTs.
#' @param project_name Name of project. Default will name the directory as "ADT-scatterplot-default".
#' @keywords plots
#' @export
#'

scatter_ADTs <- function(obj_Seurat, no_ADTs = NULL, project_name = NULL){

  if(is.null(no_ADTs)){
    ADTs <- rownames(obj_Seurat@assays$ADT@data)
  }
  else {
    remove_ADT <- paste0(no_ADTs, sep = "|")
    ADTs <- rownames(obj_Seurat@assays$ADT@data)[!grepl(remove_ADT, rownames(obj_Seurat@assays$ADT@data))]
  }

  if(is.null(project_name)){
    project_name <- "default"
  }

  DefaultAssay(obj_Seurat) <- "ADT"
  VariableFeatures(obj_Seurat) <- ADTs

  # Get all combinations of ADTs
  plot_adt <- t(combn(ADTs, 2))
  adt_data <- as.data.frame(t(as.matrix(GetAssayData(obj_Seurat[["ADT"]], slot = "data"))))

  plist <- list(theme_bw(), geom_density_2d(color = 'deepskyblue3'),
                geom_vline(xintercept = 0, linetype = 'dashed'),
                geom_hline(yintercept = 0, linetype = 'dashed')
  )

  for(i in 1:nrow(plot_adt)){
    p <- qplot(x = adt_data[,plot_adt[i,1]],
               y = adt_data[,plot_adt[i,2]],
               size = I(0.3)) +
      plist + xlab(plot_adt[i,1]) + ylab(plot_adt[i,2])
    ggsave(paste0("ADT-scatterplot-", project_name, "/", plot_adt[i,1], "vs", plot_adt[i,2], ".png"),
         plot = p, width = 3, height = 3)
  }

}
