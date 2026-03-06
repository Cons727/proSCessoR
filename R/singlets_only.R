#' Filter Singlets
#'
#' This function uses the output of demultiplex in order to filter out empty drops and doublets, keeping only singlets. It returns a Seurat object.
#' @param demux_Seurat Demultiplexed Seurat object. Output from demultiplex.
#' @keywords singlets
#' @export
#'

singlets_only <- function(demux_Seurat){
  #demux_Seurat has Idents = ident.global from previous function
  singlet_Seurat = subset(demux_Seurat, idents = "Singlet")

  singlet_Seurat
  }
