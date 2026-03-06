#' HTO Demultiplexing
#'
#' This function uses the output of bg_correction in order to demultiplex and further remove empty droplets or multiplets. It returns a Seurat object.
#' @param dsb_Seurat Seurat object. Output from bg_correction.
#' @param pool_HTO Vector of HTOs in the assay. Default all HTOs.
#' @param method Demultiplex method, HTODemux, MULTI, demuxmix. Default demuxmix.
#' @param mt_threshold Proportion of mitochondrial gene used as cutoff. Default 0.2.
#' @param make_plots Logical. If TRUE plots are generated and saved as PNG. Defaults to TRUE.
#' @param model Character specifying the type of mixture model to used in demuxmix. Either "naive", "reg", "regpos", or "auto" (default).
#' @param clusterInit Optional list of numeric vectors to manually specify the droplet to component assignment used to initialize the EM algorithm of demuxmix. The name of each list element must match a row name of hto. The length of each element must match the number of columns of hto. Only the values 1 and 2 are allowed, where 1 belongs to the negative.
#' @keywords demux
#' @export
#'

demultiplex <- function(dsb_Seurat, pool_HTO=NULL, method="demuxmix", mt_threshold=0.2, make_plots=TRUE, model="auto", clusterInit=NULL){

  print("Starting demultiplex")

  ident.id <- "hash.ID"
  ident.global <- "HTO_classification.global"
  ident.class <- "HTO_classification"

  # Flag for assay without RNA
  make_violin <- TRUE

  if(is.null(pool_HTO)){
    pool_HTO <- rownames(dsb_Seurat@assays$HTO@data)
  }

  if(method == "demuxmix"){

    if(is.null(dsb_Seurat@assays$RNA)){ # No RNA

      dsb_Seurat <- NormalizeData(dsb_Seurat, assay = "HTO", normalization.method = "CLR", verbose = FALSE,  margin = 2)
      hto <- as.matrix(GetAssayData(dsb_Seurat, slot = "counts", assay = "HTO"))
      
      if(is.null(clusterInit)){
        dmm <- demuxmix::demuxmix(hto = hto, model = "naive")
      }else{
        dmm <- demuxmix::demuxmix(hto = hto, model = "naive", clusterInit = clusterInit)
      }

      make_violin <- FALSE

    }else{

      dsb_Seurat <- subset(dsb_Seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < mt_threshold*100)
      dsb_Seurat <- NormalizeData(dsb_Seurat, assay = "HTO", normalization.method = "CLR", verbose = FALSE,  margin = 2)
      hto <- as.matrix(GetAssayData(dsb_Seurat, slot = "counts", assay = "HTO"))
      
      if(is.null(clusterInit)){
        n.genes <- dsb_Seurat@meta.data[["nFeature_RNA"]]
        dmm <- demuxmix(hto = hto, rna = n.genes, model = model)
      }else{
        dmm <- demuxmix::demuxmix(hto = hto, model = "naive", clusterInit = clusterInit)
      }
    }

    dmm <- dmmClassify(dmm)

    # Named in the same way as Seurat functions
    dsb_Seurat[[ident.global]] <- dmm %>%
      mutate(global = factor(case_when(Type == "singlet" ~ "Singlet",
                       Type == "multiplet" ~ "Multiplet",
                       Type %in% c("uncertain", "negative") ~ "Negative"))
             ) %>%
      select(global)

    dsb_Seurat[[ident.class]] <- dmm %>%
      mutate(id = case_when(Type == "singlet" ~ HTO,
                            Type %in% c("uncertain", "negative") ~ "Negative")) %>%
      select(id)
    dsb_Seurat[[ident.id]] <- dmm %>%
      mutate(id = case_when(Type == "singlet" ~ HTO,
                            Type == "multiplet" ~ "Multiplet",
                            Type %in% c("uncertain", "negative") ~ "Negative")) %>%
      select(id)
  }

  if(method == "HTODemux"){

    if(is.null(dsb_Seurat@assays$RNA)){ # No RNA
      dsb_Seurat <- subset(dsb_Seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 5000)
    } else{
      dsb_Seurat <- subset(dsb_Seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < mt_threshold*100)
    }
    dsb_Seurat <- NormalizeData(dsb_Seurat, assay = "HTO", normalization.method = "CLR", verbose = FALSE,  margin = 2)

    # Determine what is the best positive quantile for demuxing
    DemuxResults <- data.frame(Positive.Quantile = seq(0.8, 0.99, 0.01), Cell.Recovery = NA)

    pb = txtProgressBar(min = 0, max = length(DemuxResults$Positive.Quantile), initial = 0, style = 3)

    for(i in 1:length(DemuxResults$Positive.Quantile)) {
      dsb_Seurat = HTODemux(dsb_Seurat, assay = "HTO", positive.quantile =  DemuxResults$Positive.Quantile[i])
      DemuxResults$Cell.Recovery[i] = table(dsb_Seurat$hash.ID %in% pool_HTO)[[2]]
      setTxtProgressBar(pb, i)
    }
    # TODO?: Add criteria, if every HTO does not have at least 2 reads, adjust maximum cell recovery.
    PosQuant <- DemuxResults[which.max(DemuxResults$Cell.Recovery),][[1]]
    print(paste0("Using ", PosQuant, " as positive quantile"))
    dsb_Seurat <- HTODemux(dsb_Seurat, assay = "HTO", positive.quantile =  PosQuant)

    # Change Doublet to Multiplet
    dsb_Seurat@meta.data[[ident.id]] <- ifelse(dsb_Seurat[[ident.id]] == "Doublet", "Multiplet", unlist(dsb_Seurat[[ident.id]]))
    levels(dsb_Seurat@meta.data[[ident.global]])[levels(dsb_Seurat@meta.data[[ident.global]]) == "Doublet"] <- "Multiplet"
  }

  if(method == "MULTI"){

    if(is.null(dsb_Seurat@assays$RNA)){ # No RNA
      dsb_Seurat <- subset(dsb_Seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 5000)
    } else{
      dsb_Seurat <- subset(dsb_Seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < mt_threshold*100)
    }
    dsb_Seurat <- NormalizeData(dsb_Seurat, assay = "HTO", normalization.method = "CLR", verbose = FALSE,  margin = 2)

    dsb_Seurat <- MULTIseqDemux(dsb_Seurat,
                       autoThresh = TRUE,
                       maxiter = 10,
                       qrange = seq(from = 0.1, to = 0.9, by = 0.05),
                       verbose = TRUE)

    # Rename metadata
    dsb_Seurat[[ident.id]] <- dsb_Seurat[["MULTI_ID"]]
    dsb_Seurat[[ident.class]] <- dsb_Seurat[["MULTI_classification"]]

    dsb_Seurat[["MULTI_ID"]] <- NULL
    dsb_Seurat[["MULTI_classification"]] <- NULL

    dsb_Seurat[[ident.global]] <- factor(ifelse(dsb_Seurat[[ident.id]] == "Doublet", "Multiplet",
                                                                      ifelse(dsb_Seurat[[ident.id]] == "Negative", "Negative",
                                                                             "Singlet"))
                                        )
  }

  hto <- table(dsb_Seurat@meta.data$hash.ID, dsb_Seurat@meta.data$orig.ident)
  hto <- as.data.frame.matrix(hto)
  write.csv(hto, paste0("Processing_QC/", dsb_Seurat@project.name, "_dsb_HTOcounts.csv"), quote = FALSE)
  print(hto)

  if(make_plots){
    # Plot heatmap
    # if(method == "HTODemux"){
    #   Idents(dsb_Seurat) <- ident.id
    #   p <- HTOHeatmap(dsb_Seurat, assay = "HTO", raster = FALSE) + ggtitle(dsb_Seurat@project.name)
    #   ggsave(paste0("Processing_QC/",dsb_Seurat@project.name, "_dsb_HTOheatmap.png"), plot = p, width = 8, height = 8, units = "in")

    #   # Plot per batch
    #   for(i in unique(dsb_Seurat@meta.data[["orig.ident"]])){
    #     so <- subset(dsb_Seurat, subset = orig.ident == i)
    #     p <- HTOHeatmap(so, assay = "HTO", raster = FALSE) + ggtitle(i)
    #     ggsave(paste0("Processing_QC/", i, "_dsb_HTOheatmap.png"), plot = p)
    #   }

    # } else{
      # p <- HTOMultiHeatmap(object = dsb_Seurat) + ggtitle(dsb_Seurat@project.name)
      # ggsave(paste0("Processing_QC/", dsb_Seurat@project.name, "_dsb_HTOheatmap.png"), plot = p, width = 8, height = 8, units = "in")

      # Plot per batch
      for(i in unique(dsb_Seurat@meta.data[["orig.ident"]])){
        so <- subset(dsb_Seurat, subset = orig.ident == i)
        p <- HTOMultiHeatmap(object = so) + ggtitle(i)
        ggsave(paste0("Processing_QC/", i, "_dsb_HTOheatmap.png"), plot = p, width = 8, height = 8, units = "in")
      }
    #}

    # Plot Ridge plot
    Idents(dsb_Seurat) <- ident.id
    p <- RidgePlot(dsb_Seurat, assay = "HTO", features = pool_HTO, ncol = 3)
    p <- p + patchwork::plot_annotation(title = dsb_Seurat@project.name)
    ggsave(paste0("Processing_QC/", dsb_Seurat@project.name, "_dsb_HTORidgePlot.png"), plot = p, width = 15, height = 15, units = "in")

    ## Plot per batch
      for(i in unique(dsb_Seurat@meta.data[["orig.ident"]])){
        so <- subset(dsb_Seurat, subset = orig.ident == i)
        p <- RidgePlot(so, assay = "HTO", features = pool_HTO, ncol = 3)
        p <- p + patchwork::plot_annotation(title = i)
        ggsave(paste0("Processing_QC/", i, "_dsb_HTORidgePlot.png"), plot = p, width = 15, height = 15, units = "in")
      }
  
    # Plot Violin plots
    if(make_violin){
      Idents(dsb_Seurat) <- ident.global
      p <- VlnPlot(dsb_Seurat, features = "nCount_RNA", pt.size = 0, log = TRUE) + ggtitle(dsb_Seurat@project.name)
      ggsave(paste0("Processing_QC/", dsb_Seurat@project.name, "_dsb_UMIcounts.png"), plot = p)

      ## Plot per batch
      for(i in unique(dsb_Seurat@meta.data[["orig.ident"]])){
        so <- subset(dsb_Seurat, subset = orig.ident == i)
        p <- VlnPlot(so, features = "nCount_RNA", pt.size = 0, log = TRUE) + ggtitle(i)
        ggsave(paste0("Processing_QC/", i, "_dsb_UMIcounts.png"), plot = p)
      }

      Idents(dsb_Seurat) <- ident.global
      p <- VlnPlot(dsb_Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, combine = FALSE)
      new.titles <- c("Genes per Cell", "UMI per Cell", "Percent MT")
      p <- lapply(p %>% seq, function(i) {
        p[[i]] & ggtitle(new.titles[i]) & NoLegend()
      })
      p <- p[[1]] + p[[2]] + p[[3]] + patchwork::plot_annotation(title = dsb_Seurat@project.name)
      ggsave(paste0("Processing_QC/", dsb_Seurat@project.name, "_dsb_UMIcountsParsed.png"), plot = p, width = 15, height = 10, units = "in")

      ## Plot per batch
      for(i in unique(dsb_Seurat@meta.data[["orig.ident"]])){
        so <- subset(dsb_Seurat, subset = orig.ident == i)
        p <- VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, combine = FALSE)
        new.titles <- c("Genes per Cell", "UMI per Cell", "Percent MT")
        p <- lapply(p %>% seq, function(i) {
          p[[i]] & ggtitle(new.titles[i]) & NoLegend()
        })
        p <- p[[1]] + p[[2]] + p[[3]] + patchwork::plot_annotation(title = i)
        ggsave(paste0("Processing_QC/", i, "_dsb_UMIcountsParsed.png"), plot = p, width = 15, height = 10, units = "in")
      }
      
    }
    ## Plot HTO tSNE colored by HTO counts, global and HTO classification
    so <- dsb_Seurat
    DefaultAssay(so) <- "HTO"
    so <- NormalizeData(so, normalization.method = "CLR")
    VariableFeatures(so) <- rownames(so@assays$HTO)
    so <- ScaleData(so) %>%
      RunPCA(approx = FALSE) %>%
      RunTSNE(dims = 1:length(VariableFeatures(so)), perplexity = 100, check_duplicates = FALSE)

    p <- FeaturePlot(so, features = VariableFeatures(so), reduction = "tsne")
    ggsave(paste0("Processing_QC/", so@project.name, "_HTO_counts_tSNE.png"), plot = p, width = 16, height = 16)

    p <- DimPlot(so, reduction = "tsne", group.by = "HTO_classification.global")
    ggsave(paste0("Processing_QC/", so@project.name, "_HTO_global_tSNE.png"), plot = p, width = 8, height = 8)

    p <- DimPlot(so, reduction = "tsne", group.by = "HTO_classification")
    ggsave(paste0("Processing_QC/", so@project.name, "_HTO_tSNE.png"), plot = p, width = 8, height = 8)

  }

  Idents(dsb_Seurat) <- ident.global

  dsb_Seurat
}
