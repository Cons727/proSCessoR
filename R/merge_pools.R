#' Define ADT Background Threshold from Merged Pools
#'
#' This function returns the QC scatter plot and density plot to define correct thresholds for background drops from ADT assays. Returns a list of Seurat objects. It requires count/raw_feature_bc_matrix and count/sample_feature_bc_matrix from CellRanger outputs. The plots are saved to the current directory as QC_scatter_noBgThreshold.png. Define the thresholds by identifying the more dense population of background drops that is not close to 0. Normally, the background drops are located in the second peak.
#' @param batches Vector of CellRanger output file name to be processed, i.e. "BCS2024_P2_AgNeg".
#' @param mt_threshold Proportion of mitochondrial gene used as cutoff. Default 0.2.
#' @param project_name (optional) Project name used in the creation of the Seurat object.
#' @param extract_metrics Logical. Save metrics to file. Default TRUE.
#' @keywords bg
#' @export

merge_pools <- function(batches, mt_threshold=0.2, project_name=NULL, extract_metrics=TRUE){

  if(is.null(project_name)){
    project_name <- "SeuratProject"
  }

  obj_seurat_raw <- list()
  obj_seurat_filter <- list()
  stopifnot(length(batches) > 0)

  for(batch in batches){

      if(extract_metrics){
        file <- paste0(batch, "/per_sample_outs/", batch)
        file <- get_CR_file(file = "metrics_summary.csv", path = file)[1] # get file from first path
        extract_metrics(metrics_summary = file, batch = batch, to_file = TRUE)
      }

      raw_path <- paste0(batch, "/multi/count")
      raw_path <- get_CR_file(path = raw_path, file = "raw_feature_bc_matrix")[1] # get dir not .h5

      filter_path <- paste0(batch, "/per_sample_outs/", batch, "/count")
      filter_path <- get_CR_file(path = filter_path, file = "sample_filtered_feature_bc_matrix")[1] # get dir not .h5

      raw_data <- Read10X(raw_path)
      filter_data <- Read10X(filter_path)

      # Flag if data includes GEX expression, FALSE if not
      rna.exp <- any(grepl("Gene Expression", names(filter_data)))
      # Flag if data includes FBC library, FALSE if not
      fbc.exp <- any(grepl("Antibody Capture", names(filter_data)))

      if(fbc.exp){
        # Flag if data includes ADT expression, FALSE if only HTO
        adt.exp <- any(!grepl("HTO", rownames(filter_data$`Antibody Capture`)))
        # Flag if data includes HTO expression, FALSE if not
        hto.exp <- any(grepl("HTO", rownames(filter_data$`Antibody Capture`)))
      } else {
        adt.exp <- FALSE
        hto.exp <- FALSE
      }

      # RNA + ADT +/- HTO
      if(adt.exp & rna.exp & hto.exp){
        stained_cells <- colnames(filter_data$`Gene Expression`)
        background <- setdiff(colnames(raw_data$`Gene Expression`), stained_cells)

        raw_ADT <- raw_data$`Antibody Capture`[!grepl("HTO", raw_data$`Antibody Capture`@Dimnames[[1]]), ]
        raw_RNA <- raw_data$`Gene Expression`

        mtgene <- grep(pattern = "^MT-", rownames(raw_RNA), value = TRUE)

        raw_metadata <- data.frame(
          rna.size = log10(Matrix::colSums(raw_RNA)),
          prot.size = log10(Matrix::colSums(raw_ADT)),
          n.gene = Matrix::colSums(raw_RNA > 0),
          mt.prop = Matrix::colSums(raw_RNA[mtgene, ]) / Matrix::colSums(raw_RNA)
        )

        raw_metadata$drop.class = ifelse(rownames(raw_metadata) %in% stained_cells, 'cell', 'background')
        raw_metadata = raw_metadata[raw_metadata$rna.size > 0 & raw_metadata$prot.size > 0, ]

        cell_metadata <- raw_metadata[raw_metadata$drop.class == 'cell' & raw_metadata$n.gene > 200, ]
        rna.mult <- 3*mad(cell_metadata$rna.size)
        prot.mult <- 3*mad(cell_metadata$prot.size)
        rna.lower <- median(cell_metadata$rna.size) - rna.mult
        rna.upper <- median(cell_metadata$rna.size) + rna.mult
        prot.lower <- median(cell_metadata$prot.size) - prot.mult
        prot.upper <- median(cell_metadata$prot.size) + prot.mult

        qc_cells <- rownames(
          cell_metadata[cell_metadata$prot.size > prot.lower &
                          cell_metadata$prot.size < prot.upper &
                          cell_metadata$rna.size > rna.lower &
                          cell_metadata$rna.size < rna.upper &
                          cell_metadata$mt.prop < mt_threshold, ]
        )

        cell.adt.raw <- as.matrix(raw_ADT[ , qc_cells])
        cell.rna.raw <- raw_RNA[ , qc_cells]
        cell_metadata <- cell_metadata[qc_cells, ]

        obj_seurat_filter[[batch]] <- CreateSeuratObject(counts = cell.rna.raw, meta.data = cell_metadata,
                                                        project = batch, assay = "RNA", min.cells = 3, min.features = 200)
        cell.adt.raw <- cell.adt.raw[, rownames(obj_seurat_filter[[batch]]@meta.data)]

        obj_seurat_filter[[batch]][["e_barcode"]] <- Cells(obj_seurat_filter[[batch]])

        obj_seurat_filter[[batch]][["ADT"]] <- CreateAssayObject(counts = cell.adt.raw)

        # Define background
        background_drops <- rownames(raw_metadata[!rownames(raw_metadata) %in% qc_cells, ])

        background.adt.mtx <- as.matrix(raw_ADT[ , background_drops])
        background_metadata <- raw_metadata[background_drops, ]

        obj_seurat_raw[[batch]] <- CreateSeuratObject(counts = background.adt.mtx, project = batch,
                                                      assay = "ADT", min.cells = 0, min.features = 0,
                                                      meta.data = background_metadata)

        if(hto.exp){
          hto.counts <- raw_data$`Antibody Capture`[grep("HTO", raw_data$`Antibody Capture`@Dimnames[[1]]), ]

          if(is.null(nrow(hto.counts))){ # One hashtag, in some weird cases, you want to multiplex one sample
            tmp <- matrix(hto.counts[colnames(x = obj_seurat_filter[[batch]])], nrow = 1,
                          ncol = length(colnames(x = obj_seurat_filter[[batch]])))
            colnames(tmp) <- colnames(x = obj_seurat_filter[[batch]])
            rownames(tmp) <- grep("HTO", raw_data$`Antibody Capture`@Dimnames[[1]], value = TRUE)
            
            obj_seurat_filter[[batch]][["HTO"]] <- CreateAssayObject(counts = tmp)

          } else{
            obj_seurat_filter[[batch]][["HTO"]] <- CreateAssayObject(counts = hto.counts[, colnames(x = obj_seurat_filter[[batch]])])
          }
        }
      }

      # RNA +/- HTO
      if(rna.exp & !adt.exp){
        obj_seurat_filter[[batch]] <- rna_correction(raw_data = raw_data,
                                                      filter_data = filter_data,
                                                      project_name = batch,
                                                      mt_threshold = mt_threshold,
                                                      hto_exp = hto.exp)
        obj_seurat_filter[[batch]][["e_barcode"]] <- Cells(obj_seurat_filter[[batch]])
      }

      # ADT +/- HTO
      if(!rna.exp & adt.exp){
        stained_cells <- colnames(filter_data$`Antibody Capture`)

        raw_ADT <- raw_data$`Antibody Capture`[!grepl("HTO", raw_data$`Antibody Capture`@Dimnames[[1]]), ]

        raw_metadata <- data.frame(prot.size = log10(Matrix::colSums(raw_ADT)))

        raw_metadata$drop.class = ifelse(rownames(raw_metadata) %in% stained_cells, 'cell', 'background')
        raw_metadata = raw_metadata[raw_metadata$prot.size > 0, ]

        background_drops <- rownames(raw_metadata[raw_metadata$drop.class == 'background',])

        background.adt.mtx <- as.matrix(raw_ADT[ , background_drops])
        background_metadata <- raw_metadata[background_drops, ]

        obj_seurat_raw[[batch]] <- CreateSeuratObject(counts = background.adt.mtx, project = batch,
                                                      assay = "ADT", min.cells = 0, min.features = 0,
                                                      meta.data = background_metadata)

        cell_metadata <- raw_metadata[raw_metadata$drop.class == 'cell', ]
        prot.mult <- (3*mad(cell_metadata$prot.size))
        prot.lower <- median(cell_metadata$prot.size) - prot.mult
        prot.upper <- median(cell_metadata$prot.size) + prot.mult

        qc_cells <- rownames(
          cell_metadata[cell_metadata$prot.size > prot.lower &
                        cell_metadata$prot.size < prot.upper, ]
        )

        cell.adt.raw <- as.matrix(raw_ADT[ , qc_cells])
        cell_metadata <- cell_metadata[qc_cells, ]

        obj_seurat_filter[[batch]] <- CreateSeuratObject(counts = cell.rna.raw, meta.data = cell_metadata,
                                                        project = batch, assay = "RNA", min.cells = 3, min.features = 200)
        cell.adt.raw <- cell.adt.raw[, rownames(obj_seurat_filter[[batch]]@meta.data)]

        obj_seurat_filter[[batch]][["e_barcode"]] <- Cells(obj_seurat_filter[[batch]])

        obj_seurat_filter[[batch]][["ADT"]] <- CreateAssayObject(counts = cell.adt.raw)

        if(hto.exp){
          hto.counts <- raw_data$`Antibody Capture`[grep("HTO", raw_data$`Antibody Capture`@Dimnames[[1]]), ]
          obj_seurat_filter[[batch]][["HTO"]] <- CreateAssayObject(counts = hto.counts[, colnames(x = obj_seurat_filter[[batch]])])
        }
      }

      # HTO
      if(!rna.exp & !adt.exp & hto.exp){
        hto.counts <- filter_data$`Antibody Capture`[grep("HTO", filter_data$`Antibody Capture`@Dimnames[[1]]), ]
        obj_seurat_filter[[batch]][["HTO"]] <- CreateAssayObject(counts = hto.counts)
      }
    }

  # Merge ADT +/- HTO +/- RNA
  if(adt.exp){
    if(length(obj_seurat_filter) > 1){
      raw_data <- merge(obj_seurat_raw[[1]], y = obj_seurat_raw[-1], add.cell.ids = names(obj_seurat_raw))
      filter_data <- merge(obj_seurat_filter[[1]], y = obj_seurat_filter[-1], add.cell.ids = names(obj_seurat_filter), project = project_name)
    } else {
      raw_data <- obj_seurat_raw[[1]]
      filter_data <- obj_seurat_filter[[1]]
    }

    # Plot noBgThreshold RNA + ADT
    if(rna.exp){
      filter_data[["percent.mt"]] <- PercentageFeatureSet(filter_data, pattern = "^MT-")

      df <- data.frame(rna.size = c(raw_data@meta.data$rna.size, filter_data@meta.data$rna.size),
                      prot.size = c(raw_data@meta.data$prot.size, filter_data@meta.data$prot.size),
                      mt.prop = c(raw_data@meta.data$mt.prop, filter_data@meta.data$mt.prop))
      rna.lower <- min(filter_data@meta.data$rna.size)
      rna.upper <- max(filter_data@meta.data$rna.size)
      prot.lower <- min(filter_data@meta.data$prot.size)
      prot.upper <- max(filter_data@meta.data$prot.size)

      # Select background drops
      cutoffs <- c(rep(NA, 4))
      names(cutoffs) <- c("rna.low", "rna.high", "prot.low", "prot.high")

      # RNA
      model <- normalmixEM(x = filter(raw_data@meta.data, drop.class != "cell")$rna.size, k = 2)
      index.higher <- which.max(model$mu)  # Index of component with higher mean (non-zero background peak)

      cutoffs["rna.low"] <- find.cutoff(proba = 0.9999, model = model, i = index.higher)
      cutoffs["rna.high"] <- 2*model$mu[index.higher] - cutoffs["rna.low"]

      # Protein
      model <- normalmixEM(x = filter(raw_data@meta.data, drop.class != "cell")$prot.size, k = 2)
      index.higher <- which.max(model$mu)  # Index of component with higher mean (non-zero background peak)

      cutoffs["prot.low"] <- find.cutoff(proba = 0.9999, model = model, i = index.higher)
      cutoffs["prot.high"] <- 2*model$mu[index.higher] - cutoffs["prot.low"]

      # Don't take less than 0.8
      if(any(cutoffs < 0.8)){
        cutoffs[cutoffs < 0.8] <- 0.8
      }

      # Don't take more than 3
      if(any(cutoffs > 3)){
        cutoffs[cutoffs > 3] <- 3
      }

      # Make plot
      scatter <- ggplot(df, aes(x = rna.size, y = prot.size, color = mt.prop)) +
        theme_bw() +
        geom_point(size = 0.3) +
        geom_density2d(color = "skyblue") +
        scale_color_viridis_c(option = "C") +
        geom_rect(aes(xmin = rna.lower, xmax = rna.upper, ymin = prot.lower, ymax = prot.upper),
                  color = "green4", size = 0.8, fill = NA) +
        geom_rect(aes(xmin = cutoffs["rna.low"], xmax = cutoffs["rna.high"], ymin = cutoffs["prot.low"], ymax = cutoffs["prot.high"]),
                  color = "darkorange2", size = 0.8, fill = NA) +
        labs(x = "log10(RNA Counts)", y = "log10(Protein Counts)") +
        theme(legend.position = "bottom")

      densRNAx <- ggplot(df) +
        geom_density(aes(x = rna.size), fill = "#AA4AC4", color = "#AA4AC4") +
        theme_void() +
        theme(legend.position = "none") +
        labs(title = project_name)

      densProty <- ggplot(df) +
        geom_density(aes(x = prot.size), fill = "#00C1D5", color = "#00C1D5") +
        theme_void() +
        theme(legend.position = "none") +
        coord_flip()

      p <- densRNAx + patchwork::plot_spacer() + scatter + densProty +
        patchwork::plot_layout(ncol = 2, nrow = 2, widths = c(4,1), heights = c(1,4))

      ggsave(paste0("Processing_QC/", project_name, "_QC_scatter_BgThreshold.png"), plot = p, width = 5, height = 5)
    }

    # Plot ADT density
    if(!rna.exp){
      df <- data.frame(prot.size = c(raw_data@meta.data$prot.size, filter_data@meta.data$prot.size))
      prot.lower <- min(filter_data@meta.data$prot.size)
      prot.upper <- max(filter_data@meta.data$prot.size)

      cutoffs <- c(rep(NA, 2))
      names(cutoffs) <- c("prot.low", "prot.high")

      # Protein
      model <- normalmixEM(x = filter(raw_data@meta.data, drop.class != "cell")$prot.size, k = 2)
      index.higher <- which.max(model$mu)  # Index of component with higher mean (non-zero background peak)

      cutoffs["prot.low"] <- find.cutoff(proba = 0.9999, model = model, i = index.higher)
      cutoffs["prot.high"] <- 2*model$mu[index.higher] - cutoffs["prot.low"]

      ggplot(df) +
        geom_density(aes(x = prot.size), color = "#00C1D5") +
        geom_vline(aes(xintercept = prot.lower), color = "green4") +
        geom_vline(aes(xintercept = prot.upper), color = "green4") +
        geom_vline(aes(xintercept = cutoffs["prot.low"]), color = "darkorange2") +
        geom_vline(aes(xintercept = cutoffs["prot.high"]), color = "darkorange2") +
        theme(legend.position = "none") +
        labs(title = project_name)
      ggsave(paste0("Processing_QC/", project_name, "_logADT_Density.png"))
    }

    # Return seurat_list
    res <- list(raw_data = raw_data, filter_data = filter_data, mt_threshold = mt_threshold, bg_thresholds = cutoffs)
  }

  # Merge RNA +/- HTO, no plot
  if(rna.exp & !adt.exp){
    if(length(obj_seurat_filter) > 1){
      filter_data <- merge(obj_seurat_filter[[1]], y = obj_seurat_filter[-1], add.cell.ids = names(obj_seurat_filter), project = project_name)
      filter_data[["percent.mt"]] <- PercentageFeatureSet(filter_data, pattern = "^MT-")
      res <- list(filter_data = filter_data, mt_threshold = mt_threshold, bg_thresholds = NULL)
    } else{
      obj_seurat_filter[["percent.mt"]] <- PercentageFeatureSet(obj_seurat_filter, pattern = "^MT-")
      res <- list(filter_data = obj_seurat_filter, mt_threshold = mt_threshold, bg_thresholds = NULL)
    }
  }

  res

}


#' Find background threshold
#'
#' Find the background cutoff for protein normalization with dsb. Returns a numeric value.
#' @param proba Numeric value. Probability that the sample comes from the zero background component.
#' @param model List returned from normalmixEM.
#' @param i Integer. Index of the component that have the non-zero component, usually the second peak.
#' @param upper.cut Integer. Upper value to look for the root.
#' @keywords preprocessing
#' @export
#'

find.cutoff <- function(proba = 0.5, model, i = index.higher) {
  x <- model$x
  
  ## Cutoff such that Pr[drawn from non-zero background component] == proba
  f <- function(x) {
    proba - (model$lambda[i]*pnorm(x, model$mu[i], model$sigma[i]) /
               (model$lambda[1]*pnorm(x, model$mu[1], model$sigma[1]) + model$lambda[2]*pnorm(x, model$mu[2], model$sigma[2])))
    }
  root <- tryCatch({
    uniroot(f=f, lower = 0, upper = 5)$root
  },
  error = function(e) {
    message('An error occurred when estimating cutoff, setting lower threshold = 1')
    return(1)
  })
  return(root)
}
