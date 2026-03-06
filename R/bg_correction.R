#' RNA+ADT Background Noise Correction/Removal
#'
#' This function removes background noise from ADT assays using DSB package and filters droplets based on RNA counts and mitochondrial proportion. It requires count/raw_feature_bc_matrix and count/sample_feature_bc_matrix from CellRanger outputs. Returns a Seurat object with denoised ADT counts in data slot of ADT assay.
#' If isotype controls were not included in the experiment, function results in not defining the technical component of each cell’s protein library. The values of the normalized matrix returned are the number of standard deviations above the expected ambient noise captured by empty droplets. It is strongly recommend to use isotype controls. However if these are not available, the background mean for each cell inferred via a per-cell Gaussian mixture model can theoretically be used alone to define the cell’s technical component, but this assumes the background mean has no expected biological variation. Isotype controls anchor the component of the background mean associated with noise.
#' If the normalize counts have a large range, select the correct background threshold by running bg_thresholds function and manually defining the values.
#' @param raw_data CellRanger output raw_feature_bc_matrix.
#' @param filter_data CellRanger output sample_feature_bc_matrix.
#' @param isotype Name of the isotype used as ADT control. Default to HSA.
#' @param denoise_bg Denoising with background mean only. Default to TRUE.
#' @param pool_HTO Vector of HTOs in the assay to demultiplex. Default all HTOs.
#' @param project_name (optional) Project name used in the creation of the Seurat object.
#' @param mt_threshold Proportion of mitochondrial gene used as cutoff. Default 0.2.
#' @param prot.bg.threshold Vector of protein size thresholds to define drops as background. Default c(1.5, 3).
#' @param rna.bg.threshold Vector of RNA size thresholds to define drops as background. Default c(0, 2.5).
#' @keywords bg
#' @export
#'

bg_correction <- function(raw_data, filter_data,
                          isotype="HSA", denoise_bg=TRUE,
                          pool_HTO=NULL,
                          project_name=NULL,
                          mt_threshold=0.2, prot.bg.threshold=c(1.5, 3), rna.bg.threshold=c(0, 2.5)){

  print("Starting background removal")

  if(is.null(project_name)){
    project_name <- "SeuratProject"
  }

  stained_cells <- colnames(filter_data$`Gene Expression`)
  # Split assays and remove HTO from ADT
  raw_ADT <- raw_data$`Antibody Capture`[!grepl("HTO", raw_data$`Antibody Capture`@Dimnames[[1]]), ]
  raw_RNA <- raw_data$`Gene Expression`

  ## Identify "empty" drops
  # Mitochondrial reads
  mtgene <- grep(pattern = "^MT-", rownames(raw_RNA), value = TRUE)

  raw_metadata <- data.frame(
    rna.size = log10(Matrix::colSums(raw_RNA)),
    prot.size = log10(Matrix::colSums(raw_ADT)),
    n.gene = Matrix::colSums(raw_RNA > 0),
    mt.prop = Matrix::colSums(raw_RNA[mtgene, ]) / Matrix::colSums(raw_RNA)
  )

  # Add an indicator for which barcodes CellRanger called as cells
  raw_metadata$drop.class = ifelse(rownames(raw_metadata) %in% stained_cells, 'cell', 'background')

  # Remove barcodes with no evidence of capture in the experiment
  raw_metadata = raw_metadata[raw_metadata$rna.size > 0 & raw_metadata$prot.size > 0, ]

  #write.csv(raw_metadata, "raw_metadata.csv", row.names = FALSE)

  # Define upper and lower cutoffs to identify background drops and remove them from ADT assay
  background_drops <- rownames(
    raw_metadata[raw_metadata$prot.size > prot.bg.threshold[1] &
                   raw_metadata$prot.size < prot.bg.threshold[2] &
                   raw_metadata$rna.size > rna.bg.threshold[1] &
                   raw_metadata$rna.size < rna.bg.threshold[2], ]
  )

  background.adt.mtx <- as.matrix(raw_ADT[ , background_drops])

  # Calculate statistical thresholds for droplet filtering.
  cell_metadata <- raw_metadata[raw_metadata$drop.class == 'cell', ]

  # Filter drops with +/- 3 median absolute deviations from the median library size
  rna.mult <- (3*mad(cell_metadata$rna.size))
  prot.mult <- (3*mad(cell_metadata$prot.size))
  rna.lower <- median(cell_metadata$rna.size) - rna.mult
  rna.upper <- median(cell_metadata$rna.size) + rna.mult
  prot.lower <- median(cell_metadata$prot.size) - prot.mult
  prot.upper <- median(cell_metadata$prot.size) + prot.mult

  cat(paste0("Threshold for filtering drops as cells\nRNA lower: ", rna.lower,
               "\nRNA upper: ", rna.upper,
               "\nADT lower: ", prot.lower,
               "\nADT upper: ", prot.upper), "\n")

  # Plot QC_scatter
  p2 <- ggplot(raw_metadata, aes(x = rna.size, y = prot.size , color = mt.prop)) +
    theme_bw() +
    geom_point(size = 0.3) +
    scale_color_viridis_c(option = "C") +
    geom_rect(aes(xmin = rna.bg.threshold[1], xmax = rna.bg.threshold[2], ymin = prot.bg.threshold[1], ymax = prot.bg.threshold[2]),
              color = "skyblue", size = 0.8, fill = NA) +
    geom_rect(aes(xmin = rna.lower, xmax = rna.upper, ymin = prot.lower, ymax = prot.upper),
              color = "darkorange2", size = 0.8, fill = NA) +
    labs(title = project_name, x = "log10(RNA Counts)", y = "log10(Protein Counts)")
  ggsave("QC_scatter_bgThreshold.png", plot = p2)

  # Filter rows based on droplet quality control metrics
  qc_cells <- rownames(
    cell_metadata[cell_metadata$prot.size > prot.lower &
                    cell_metadata$prot.size < prot.upper &
                    cell_metadata$rna.size > rna.lower &
                    cell_metadata$rna.size < rna.upper &
                    cell_metadata$mt.prop < mt_threshold, ]
  )

  # Keep only QC cells in metadata, ADT and RNA
  cell.adt.raw <- as.matrix(raw_ADT[ , qc_cells])
  cell.rna.raw <- raw_RNA[ , qc_cells]
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
  stopifnot(isTRUE(all.equal(rownames(cell_metadata), colnames(cell.rna.raw))))

  dsb_Seurat <- CreateSeuratObject(counts = cell.rna.raw, meta.data = cell_metadata,
                                   project = project_name, assay = "RNA", min.cells = 3, min.features = 200)

  # Filter the cells_dsb_norm data to include only those cells that were included in the Seurat object creation
  cells_dsb_norm <- cells_dsb_norm[, row.names(dsb_Seurat@meta.data)]

  dsb_Seurat[["ADT"]] <- CreateAssayObject(data = cells_dsb_norm)

  hto.counts <- raw_data$`Antibody Capture`[grep("HTO", raw_data$`Antibody Capture`@Dimnames[[1]]), ]

  if(!is.null(pool_HTO)){
    hto.counts <- hto.counts[pool_HTO, ]
  }

  hto.counts <- hto.counts[ , colnames(cells_dsb_norm)]

  dsb_Seurat[["HTO"]] <- CreateAssayObject(counts = hto.counts[, colnames(x = dsb_Seurat)])

  dsb_Seurat[["percent.mt"]] <- PercentageFeatureSet(dsb_Seurat, pattern = "^MT-")

  dsb_Seurat
}
