#' Define ADT Background Threshold
#'
#' This function returns the QC scatter plot and density plot to define correct thresholds for background drops from ADT assays. It requires count/raw_feature_bc_matrix and count/sample_feature_bc_matrix from CellRanger outputs. The plots are saved to the current directory as QC_scatter_noBgThreshold.png and QC_density.png. Define the thresholds by identifying the more dense population of background drops that is not close to 0. The density plot helps identify the RNA threshold. Normally, the background drops are located in the second pick.
#' @param raw_data Path to CellRanger output raw_feature_bc_matrix.
#' @param filter_data Path to CellRanger output sample_feature_bc_matrix.
#' @param mt_threshold Proportion of mitochondrial gene used as cutoff. Default 0.2.
#' @param project_name (optional) Project name used in the creation of the Seurat object.
#' @keywords bg
#' @export
#'

bg_thresholds <- function(raw_data, filter_data, mt_threshold=0.2, project_name=NULL){

  if(is.null(project_name)){
    project_name <- "SeuratProject"
  }

  raw_data <- Read10X(raw_data)
  filter_data <- Read10X(filter_data)

  stained_cells <- colnames(filter_data$`Gene Expression`)
  background <- setdiff(colnames(raw_data$`Gene Expression`), stained_cells)
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

  # Calculate statistical thresholds for droplet filtering.
  cell_metadata <- raw_metadata[raw_metadata$drop.class == 'cell', ]

  # Filter drops with +/- 3 median absolute deviations from the median library size
  rna.mult <- (3*mad(cell_metadata$rna.size))
  prot.mult <- (3*mad(cell_metadata$prot.size))
  rna.lower <- median(cell_metadata$rna.size) - rna.mult
  rna.upper <- median(cell_metadata$rna.size) + rna.mult
  prot.lower <- median(cell_metadata$prot.size) - prot.mult
  prot.upper <- median(cell_metadata$prot.size) + prot.mult

  # Plot QC_scatter
  scatter <- ggplot(raw_metadata, aes(x = rna.size, y = prot.size, color = mt.prop)) +
    theme_bw() +
    geom_point(size = 0.3) +
    scale_color_viridis_c(option = "C") +
    geom_rect(aes(xmin = rna.lower, xmax = rna.upper, ymin = prot.lower, ymax = prot.upper),
              color = "darkorange2", size = 0.8, fill = NA) +
    labs(x = "log10(RNA Counts)", y = "log10(Protein Counts)") +
    theme(legend.position = "bottom")

  densRNAx <- ggplot(raw_metadata) +
    geom_density(aes(x = rna.size), fill = "#AA4AC4", color = "#AA4AC4") +
    theme_void() +
    theme(legend.position = "none") +
    labs(title = project_name)

  densProty <- ggplot(raw_metadata) +
    geom_density(aes(x = prot.size), fill = "#00C1D5", color = "#00C1D5") +
    theme_void() +
    theme(legend.position = "none") +
    coord_flip()

  p <- densRNAx + plot_spacer() + scatter + densProty +
    plot_layout(ncol = 2, nrow = 2, widths = c(4,1), heights = c(1,4))

  ggsave(paste0(project_name,"_QC_scatter_noBgThreshold.png"), plot = p)

}
