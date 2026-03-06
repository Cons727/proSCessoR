#' Perform PCA
#'
#' This function normalizes, finds variable features and scales data using SCTransform or NormalizeData from Seurat, and runs PCA. It returns a Seurat object.
#' @param obj_Seurat Seurat object.
#' @param assay Assay to perform the analysis on.
#' @param norm_method Normalization method ("sct" or "lognorm"). Defaults to "sct".
#' @param vars_regress Variables to regress with SCTransform; variable from obj_Seurat metadata. Defaults to mitochondrial genes.
#' @param pca_name Name for PCA results stored in the Seurat object. Default: If RNA assay, pca_name="pca".
#' @keywords pca
#' @export
#'

seurat_PCA <- function(obj_Seurat, assay, norm_method="sct", vars_regress="percent.mt", pca_name=NULL){
  if(is.null(pca_name)){
    if(assay == "RNA"){
      pca_name <- "pca"
    }
    else {
      pca_name <- paste0(tolower(assay), ".pca")
    }
  }

  if(norm_method == "sct"){
    obj_Seurat <- SCTransform(obj_Seurat, assay = assay, vars.to.regress = vars_regress, verbose = FALSE)
    obj_Seurat <- RunPCA(obj_Seurat, assay = assay, features = VariableFeatures(obj_Seurat), reduction.name = pca_name, verbose = FALSE)
  }

  if(norm_method == "lognorm"){

    obj_Seurat <- NormalizeData(obj_Seurat, assay = assay, verbose = FALSE)
    obj_Seurat <- FindVariableFeatures(obj_Seurat, assay = assay, selection.method = "mean.var.plot", verbose = FALSE)
    obj_Seurat <- ScaleData(obj_Seurat, assay = assay, features = VariableFeatures(obj_Seurat), verbose = FALSE)
    obj_Seurat <- RunPCA(obj_Seurat, assay = assay, features = VariableFeatures(obj_Seurat), reduction.name = pca_name, verbose = FALSE)
  }
  else {
    stop("Invalid norm_method. Only sct or lognorm")
  }

  obj_Seurat

}
