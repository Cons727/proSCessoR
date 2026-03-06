#' RNA+ADT Background Noise Correction/Removal from Merged Pools
#'
#' This function removes background noise from ADT assays using DSB package and filters droplets based on RNA counts and mitochondrial proportion. It requires Seurat object list from merge_pools. Returns a Seurat object with denoised ADT counts in data slot of ADT assay.
#' If isotype controls were not included in the experiment, function results in not defining the technical component of each cell’s protein library. The values of the normalized matrix returned are the number of standard deviations above the expected ambient noise captured by empty droplets. It is strongly recommend to use isotype controls. However if these are not available, the background mean for each cell inferred via a per-cell Gaussian mixture model can theoretically be used alone to define the cell’s technical component, but this assumes the background mean has no expected biological variation. Isotype controls anchor the component of the background mean associated with noise.
#' If the normalize counts have a large range, select the correct background threshold by running bg_thresholds function and manually defining the values.
#' @param seurat_list Seurat list output from merge_pools function.
#' @param isotype Name of the isotype used as ADT control. Default to HSA.
#' @param denoise_bg Denoising with background mean only. Default to TRUE.
#' @param pool_HTO Vector of HTOs in the assay to demultiplex. Default all HTOs.
#' @param project_name (optional) Project name used in the creation of the Seurat object.
#' @param prot_bg_threshold Vector of protein size thresholds to define drops as background. Default c(1.5, 3).
#' @param rna_bg_threshold Vector of RNA size thresholds to define drops as background. Default c(1, 2.5).
#' @keywords bg
#' @export
#'

bg_correction_merged <- function(seurat_list,
                          isotype=FALSE, denoise_bg=TRUE,
                          pool_HTO=NULL,
                          project_name=NULL,
                          prot_bg_threshold=c(1.5, 3), rna_bg_threshold=c(1, 2.5)){

  if(is.null(project_name)){
    project_name <- "SeuratProject"
  }

  raw_data <- seurat_list$raw_data
  filter_data <- seurat_list$filter_data
  mt_threshold <- seurat_list$mt_threshold

  # Plot BgThreshold
  # df <- data.frame(rna.size = c(raw_data@meta.data$rna.size, filter_data@meta.data$rna.size),
  #                  prot.size = c(raw_data@meta.data$prot.size, filter_data@meta.data$prot.size),
  #                  mt.prop = c(raw_data@meta.data$mt.prop, filter_data@meta.data$mt.prop))
  # rna.lower <- min(filter_data@meta.data$rna.size)
  # rna.upper <- max(filter_data@meta.data$rna.size)
  # prot.lower <- min(filter_data@meta.data$prot.size)
  # prot.upper <- max(filter_data@meta.data$prot.size)

  # scatter <- ggplot(df, aes(x = rna.size, y = prot.size, color = mt.prop)) +
  #   theme_bw() +
  #   geom_point(size = 0.3) +
  #   scale_color_viridis_c(option = "C") +
  #   geom_density2d(color = "skyblue") +
  #   geom_rect(aes(xmin = rna_bg_threshold[1], xmax = rna_bg_threshold[2], ymin = prot_bg_threshold[1], ymax = prot_bg_threshold[2]),
  #             color = "orangered", size = 0.8, fill = NA) +
  #   geom_rect(aes(xmin = rna.lower, xmax = rna.upper, ymin = prot.lower, ymax = prot.upper),
  #             color = "green4", size = 0.8, fill = NA) +
  #   labs(x = "log10(RNA Counts)", y = "log10(Protein Counts)") +
  #   theme(legend.position = "bottom")

  # densRNAx <- ggplot(df) +
  #   geom_density(aes(x = rna.size), fill = "#AA4AC4", color = "#AA4AC4") +
  #   theme_void() +
  #   theme(legend.position = "none") +
  #   labs(title = project_name)

  # densProty <- ggplot(df) +
  #   geom_density(aes(x = prot.size), fill = "#00C1D5", color = "#00C1D5") +
  #   theme_void() +
  #   theme(legend.position = "none") +
  #   coord_flip()

  # p <- densRNAx + patchwork::plot_spacer() + scatter + densProty +
  #   patchwork::plot_layout(ncol = 2, nrow = 2, widths = c(4,1), heights = c(1,4))

  # ggsave(paste0("Processing_QC/", project_name,"_QC_scatter_BgThreshold.png"), plot = p, width = 7, height = 7)


  print("Starting background removal")

  # Select bg and cells
  rna.exp <- any(grepl("RNA", names(seurat_list[["filter_data"]]@assays))) # RNA flag, TRUE if RNA is present

  expr <- FetchData(raw_data, vars = c("prot.size", "rna.size"))

  if(rna.exp){
    raw_data <- raw_data[, which(x = (expr$prot.size > prot_bg_threshold[1] & expr$prot.size < prot_bg_threshold[2]) &
                        (expr$rna.size > rna_bg_threshold[1] & expr$rna.size < rna_bg_threshold[2])
              )]
  } else {
    raw_data <- raw_data[, which(x = (expr$prot.size > prot_bg_threshold[1] & expr$prot.size < prot_bg_threshold[2]))]
  }

  background.adt.mtx <- GetAssayData(raw_data, slot = "counts", assay = "ADT")
  cell.adt.raw <- GetAssayData(filter_data, slot = "counts", assay = "ADT")

  ## DSB
  if (!isFALSE(isotype)) {
    cells_dsb_norm <- DSBNormalizeProtein(
      cell_protein_matrix = cell.adt.raw,
      empty_drop_matrix = background.adt.mtx,
      denoise.counts = TRUE,
      use.isotype.control = TRUE,
      isotype.control.name.vec = isotype,
      quantile.clipping = TRUE,
      quantile.clip = c(0.001, 0.999))
  }
  else {
    # Denoise with background mean only
    if(denoise_bg){
      cells_dsb_norm <- DSBNormalizeProtein(
        cell_protein_matrix = cell.adt.raw,
        empty_drop_matrix = background.adt.mtx,
        denoise.counts = TRUE,
        use.isotype.control = FALSE,
        quantile.clipping = TRUE,
        quantile.clip = c(0.001, 0.999))
    }
    # Do not denoise each cell's technical component
    else {
    cells_dsb_norm <- DSBNormalizeProtein(
      cell_protein_matrix = cell.adt.raw,
      empty_drop_matrix = background.adt.mtx,
      denoise.counts = FALSE,
      quantile.clipping = TRUE,
      quantile.clip = c(0.001, 0.999))
    }
  }

  filter_data <- SetAssayData(filter_data, assay = "ADT", slot = "data", new.data = as.sparse(cells_dsb_norm))

  filter_data
}
