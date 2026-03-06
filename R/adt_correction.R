#' ADT Background Noise Correction/Removal
#'
#' This function removes background noise from ADT assays using DSB package; no RNA filtering. It requires count/raw_feature_bc_matrix and count/sample_feature_bc_matrix from CellRanger outputs. Returns a Seurat object with denoised ADT counts in data slot of ADT assay.
#' If isotype controls were not included in the experiment, function results in not defining the technical component of each cell’s protein library. The values of the normalized matrix returned are the number of standard deviations above the expected ambient noise captured by empty droplets. It is strongly recommend to use isotype controls. However if these are not available, the background mean for each cell inferred via a per-cell Gaussian mixture model can theoretically be used alone to define the cell’s technical component, but this assumes the background mean has no expected biological variation. Isotype controls anchor the component of the background mean associated with noise.
#' If the normalize counts have a large range, select the correct background threshold by running bg_thresholds function and manually defining the values.
#' @param raw_data CellRanger output raw_feature_bc_matrix.
#' @param filter_data CellRanger output sample_feature_bc_matrix.
#' @param isotype Name of the isotype used as ADT control. Default to HSA.
#' @param denoise_bg Denoising with background mean only. Default to TRUE.
#' @param pool_HTO Vector of HTOs in the assay to demultiplex. Default all HTOs.
#' @param project_name (optional) Project name used in the creation of the Seurat object.
#' @param prot.bg.threshold Vector of protein size thresholds to define drops as background. Default c(1.5, 3).
#' @keywords bg
#' @export
#'

adt_correction <- function(raw_data, filter_data,
                          isotype="HSA", denoise_bg=TRUE,
                          pool_HTO=NULL,
                          project_name=NULL,
                          prot.bg.threshold=c(1.5, 3)){

  print("Starting background removal")

  if(is.null(project_name)){
    project_name <- "SeuratProject"
  }

  stained_cells <- colnames(filter_data$`Antibody Capture`)
  raw_ADT <- raw_data$`Antibody Capture`[!grepl("HTO", raw_data$`Antibody Capture`@Dimnames[[1]]), ]

  ## Identify "empty" drops
  raw_metadata <- data.frame(prot.size = log10(Matrix::colSums(raw_ADT))
                             )

  # Add an indicator for which barcodes CellRanger called as cells
  raw_metadata$drop.class = ifelse(rownames(raw_metadata) %in% stained_cells, 'cell', 'background')

  # Density ADT size
  ggplot(raw_metadata) +
    geom_density(aes(x = prot.size), color = "#00C1D5") +
    theme(legend.position = "none") +
    labs(title = project_name)
  ggsave(paste0(project_name, "_logADT_Density.png"))

  # Remove barcodes with no evidence of capture in the experiment
  raw_metadata = raw_metadata[raw_metadata$prot.size > 0, ]

  # Define upper and lower cutoffs to identify background drops and remove them from ADT assay
  background_drops <- rownames(
    raw_metadata[raw_metadata$prot.size > prot.bg.threshold[1] &
                   raw_metadata$prot.size < prot.bg.threshold[2], ]
  )

  background.adt.mtx <- as.matrix(raw_ADT[ , background_drops])

  # Calculate statistical thresholds for droplet filtering.
  cell_metadata <- raw_metadata[raw_metadata$drop.class == 'cell', ]

  # Filter drops with +/- 3 median absolute deviations from the median library size
  prot.mult <- (3*mad(cell_metadata$prot.size))
  prot.lower <- median(cell_metadata$prot.size) - prot.mult
  prot.upper <- median(cell_metadata$prot.size) + prot.mult

  cat(paste0("Threshold for filtering drops as cells\nADT lower: ", prot.lower,
             "\nADT upper: ", prot.upper), "\n")

  # Filter rows based on droplet quality control metrics
  qc_cells <- rownames(
    cell_metadata[cell_metadata$prot.size > prot.lower &
                    cell_metadata$prot.size < prot.upper &
                    cell_metadata$rna.size > rna.lower &
                    cell_metadata$rna.size < rna.upper &
                    cell_metadata$mt.prop < mt_threshold, ]
  )

  # Keep only QC cells in metadata, ADT
  cell.adt.raw <- as.matrix(raw_ADT[ , qc_cells])
  cell_metadata <- cell_metadata[qc_cells, ]

  ## DSB
  if (isotype != FALSE) {
    cells_dsb_norm <- DSBNormalizeProtein(
      cell_protein_matrix = cell.adt.raw,
      empty_drop_matrix = background.adt.mtx,
      denoise.counts = TRUE,
      use.isotype.control = TRUE,
      isotype.control.name.vec = isotype)
  }
  else {
    # Denoise with background mean only
    if(denoise_bg){
      cells_dsb_norm <- DSBNormalizeProtein(
        cell_protein_matrix = cell.adt.raw,
        empty_drop_matrix = background.adt.mtx,
        denoise.counts = TRUE,
        use.isotype.control = FALSE)
    }
    # Do not denoise each cell's technical component
    else {
      cells_dsb_norm <- DSBNormalizeProtein(
        cell_protein_matrix = cell.adt.raw,
        empty_drop_matrix = background.adt.mtx,
        denoise.counts = FALSE)
    }
  }

  # Check the cell barcodes
  stopifnot(isTRUE(all.equal(rownames(cell_metadata), colnames(cell.adt.raw))))

  dsb_Seurat[["ADT"]] <- CreateSeuratObject(data = cells_dsb_norm, meta.data = cell_metadata,
                                            project = project_name, assay = "ADT", min.cells = 0, min.features = 0)

  hto.counts <- raw_data$`Antibody Capture`[grep("HTO", raw_data$`Antibody Capture`@Dimnames[[1]]), ]

  if(!is.null(pool_HTO)){
    hto.counts <- hto.counts[pool_HTO, ]
  }

  hto.counts <- hto.counts[ , colnames(cells_dsb_norm)]

  dsb_Seurat[["HTO"]] <- CreateAssayObject(counts = hto.counts[, colnames(x = dsb_Seurat)])

  dsb_Seurat
}
