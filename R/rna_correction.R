#' RNA Background Noise Correction/Removal
#'
#' This filters droplets based on RNA counts and mitochondrial proportion. It requires count/raw_feature_bc_matrix and count/sample_feature_bc_matrix from CellRanger outputs. Returns a Seurat object with denoised ADT counts in data slot of ADT assay.
#' @param raw_data CellRanger output raw_feature_bc_matrix.
#' @param filter_data CellRanger output sample_feature_bc_matrix.
#' @param pool_HTO Vector of HTOs in the assay to demultiplex. Default all HTOs.
#' @param project_name (optional) Project name used in the creation of the Seurat object.
#' @param mt_threshold Proportion of mitochondrial gene used as cutoff. Default 0.2.
#' @param hto_exp Logical. False if HTO expression is not included. Default to TRUE. 
#' @keywords bg
#' @export
#'

rna_correction <- function(raw_data, filter_data,
                          pool_HTO=NULL,
                          project_name=NULL,
                          mt_threshold=0.2,
                          hto_exp = TRUE){

  print("Starting background removal")

  if(is.null(project_name)){
    project_name <- "SeuratProject"
  }

  stained_cells <- colnames(filter_data$`Gene Expression`)
  raw_RNA <- raw_data$`Gene Expression`

  ## Identify "empty" drops
  # Mitochondrial reads
  mtgene <- grep(pattern = "^MT-", rownames(raw_RNA), value = TRUE)

  raw_metadata <- data.frame(
    rna.size = log10(Matrix::colSums(raw_RNA)),
    n.gene = Matrix::colSums(raw_RNA > 0),
    mt.prop = Matrix::colSums(raw_RNA[mtgene, ]) / Matrix::colSums(raw_RNA)
  )

  # Add an indicator for which barcodes CellRanger called as cells
  raw_metadata$drop.class = ifelse(rownames(raw_metadata) %in% stained_cells, 'cell', 'background')

  # # Density RNA size
  ggplot(raw_metadata) +
    geom_density(aes(x = rna.size), color = "#AA4AC4") +
    theme(legend.position = "none") +
    labs(title = project_name)
  ggsave(paste0("Processing_QC/", project_name, "_logRNA_Density.png"))

  # Remove barcodes with no evidence of capture in the experiment
  raw_metadata <- raw_metadata[raw_metadata$rna.size > 0, ]

  # Calculate statistical thresholds for droplet filtering.
  cell_metadata <- raw_metadata[raw_metadata$drop.class == 'cell', ]

  # Filter drops with +/- 3 median absolute deviations from the median library size
  rna.mult <- (3*mad(cell_metadata$rna.size))
  rna.lower <- median(cell_metadata$rna.size) - rna.mult
  rna.upper <- median(cell_metadata$rna.size) + rna.mult

  cat(paste0("Threshold for filtering drops as cells\nRNA lower: ", rna.lower,
             "\nRNA upper: ", rna.upper, "\n"))

  # Filter rows based on droplet quality control metrics
  qc_cells <- rownames(
    cell_metadata[cell_metadata$rna.size > rna.lower &
                    cell_metadata$rna.size < rna.upper &
                    cell_metadata$mt.prop < mt_threshold, ]
  )

  # Keep only QC cells in metadata, RNA
  cell.rna.raw <- raw_RNA[ , qc_cells]
  cell_metadata <- cell_metadata[qc_cells, ]

  # Check the cell barcodes
  stopifnot(isTRUE(all.equal(rownames(cell_metadata), colnames(cell.rna.raw))))

  dsb_Seurat <- CreateSeuratObject(counts = cell.rna.raw, meta.data = cell_metadata,
                                   project = project_name, assay = "RNA", min.cells = 3, min.features = 200)

  if(hto_exp){
    hto.counts <- raw_data$`Antibody Capture`

    if(!is.null(pool_HTO)){
      hto.counts <- hto.counts[pool_HTO, ]
    }
    # If only one hashtag
    if(dim(hto.counts)[1] == 1){
      hto.counts <- matrix(hto.counts[, colnames(x = dsb_Seurat)], 1, dimnames = list(rownames(hto.counts), colnames(x = dsb_Seurat)))
      dsb_Seurat[["HTO"]] <- CreateAssayObject(counts = hto.counts)

    }else{
      dsb_Seurat[["HTO"]] <- CreateAssayObject(counts = hto.counts[, colnames(x = dsb_Seurat)])
    }
  }

  dsb_Seurat[["percent.mt"]] <- PercentageFeatureSet(dsb_Seurat, pattern = "^MT-")

  dsb_Seurat
}
