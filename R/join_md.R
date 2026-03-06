#' Join Metadata to Seurat Object
#'
#' Join metadata file (must include HTO and Batch_Pool) and remove HTO from file. Saves results to dsb_sample_classifications.tsv to the current directory. Also saves and return clean singlets Seurat object.
#' @param md_path Path to metadata.csv.
#' @param demux_seurat Output from demultiplex function.
#' @param SAC_HTO HTO which includes SAC as control. Default NULL.
#' @param print_res Logical. Make HTO (PTID_visit_sort-pop) heatmaps per batch and table of cell counts per HTO. Default TRUE.
#' @keywords preprocessing
#' @export
#'

join_md <- function(md_path, demux_seurat, SAC_HTO=NULL, print_res = TRUE){
  md <- read.csv(md_path)
  # Remove possible leading and trailing white spaces
  md <- as.data.frame(apply(md, 2, str_trim))

  # status <- "HTO_classification.global"
  # hto_asgmt <- "hash.ID"

  if(is.null(SAC_HTO)){
    SAC_HTO = ""
  }


  # If only one batch, add Batch_Pool to barcode
  if(!grepl(demux_seurat@meta.data$orig.ident[1], rownames(demux_seurat@meta.data)[1], fixed = TRUE)){
    # Keep original barcodes in e_barcode
    demux_seurat[["e_barcode"]]<- rownames(demux_seurat@meta.data)
    demux_seurat <- RenameCells(demux_seurat, add.cell.id = demux_seurat@meta.data$orig.ident[1])
  }

  demux_seurat[["barcode"]] <- rownames(demux_seurat@meta.data)

  demux_seurat@meta.data <- demux_seurat@meta.data %>%
    as.data.frame() %>%
    mutate(Status = HTO_classification.global,
          HTO_Assignment = hash.ID) %>%
    dplyr::rename(Batch_Pool = orig.ident) %>%
    left_join(., md, by = c(HTO_Assignment = "HTO", Batch_Pool = "Batch_Pool"))
  rownames(demux_seurat@meta.data) <- demux_seurat@meta.data$barcode

  # For enclone (proSCessoR v1)
  demux_seurat@meta.data %>%
    select(e_barcode, HTO_classification.global, HTO_Assignment, Batch_Pool) %>%
    dplyr::rename(barcode = e_barcode, Status = HTO_classification.global) %>%
    filter(HTO_Assignment != "Negative" , !HTO_Assignment %in% SAC_HTO) %>%
    write.table(paste0("Processing_QC/", demux_seurat@project.name, "_dsb_sample_classifications.tsv"), quote = F, sep = '\t', row.names = F, col.names = T)

  # For IgBLAST (proSCessoR v2)
  singlets.barcodes <- demux_seurat@meta.data %>%
    select(e_barcode, HTO_classification.global, HTO_Assignment, Batch_Pool) %>%
    filter(!HTO_Assignment %in% SAC_HTO, HTO_classification.global == "Singlet")
  for(batch in batches){
    filter(singlets.barcodes, Batch_Pool == batch) %>% select(e_barcode) %>%
    write.table(paste0("Processing_QC/", batch, "_singlets_barcodes.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
  }

  if(print_res){
    PTID_small <- paste0("x", stringr::str_sub(demux_seurat@meta.data$PTID, start = -4))
    demux_seurat[["HTO_name"]] <- case_when(demux_seurat@meta.data$HTO_Assignment %in% SAC_HTO ~ demux_seurat@meta.data$SortPop_Small,
                                            demux_seurat@meta.data$HTO_Assignment == "Negative" ~ "Negative",
                                            demux_seurat@meta.data$HTO_Assignment == "Multiplet" ~ "Multiplet",
                                            .default = paste0(PTID_small, "_v", demux_seurat@meta.data$Visit, "_", demux_seurat@meta.data$SortPop_Small))

    # Plot per batch
    for(i in unique(demux_seurat@meta.data[["Batch_Pool"]])){
      so <- subset(demux_seurat, subset = Batch_Pool == i)

      # Change assays$HTO rownames to HTO_name, save to new assay

      hto_names <- unique(so@meta.data[["HTO_name"]])
      names(hto_names) <- unique(so@meta.data[["HTO_Assignment"]])
      hto_names <- grep("Negative|Multiplet", hto_names, value = TRUE, invert = TRUE)

      hto_names <- hto_names[match(rownames(so@assays$HTO@data), names(hto_names))]
      # Remove HTOs not present in the data (HTOs with 0 counts)
      hto_names <- na.omit(hto_names)
      so[["PTID_HTO"]] <- CreateAssayObject(data = GetAssayData(so, assay = "HTO", slot = "data")[names(hto_names),])
      rownames(so@assays$PTID_HTO@data) <- hto_names
      
      Idents(so) <- "HTO_name"
      p <- HTOMultiHeatmap(object = so, classification = "HTO_name", assay = "PTID_HTO") + ggtitle(i)
      ggsave(paste0("Processing_QC/", i, "_PTID-SPheatmap.png"), plot = p, width = 10, height = 8, units = "in")
    }

    # Save table
    #hto <- as.data.frame(table(demux_seurat@meta.data$HTO_name))
    #hto <- data.frame(Number_of_singlets = hto$Freq, row.names = hto$Var1)

    hto <- demux_seurat@meta.data %>%
      group_by(Batch_Pool, HTO_name) %>%
      summarise(Number_of_singlets = n())

    write.csv(hto, paste0("Processing_QC/", demux_seurat@project.name, "_PTID-SPcounts.csv"), quote = FALSE)
    print(hto)
  }

  # Remove Doublets
  Idents(demux_seurat) <- "Status"
  singlets_seurat <- subset(demux_seurat, idents = "Multiplet", invert = TRUE)
  invisible(ifelse(!dir.exists(file.path(".", "Seurat_objects")), dir.create(file.path(".", "Seurat_objects")), FALSE))
  saveRDS(singlets_seurat, paste0("Seurat_objects/", demux_seurat@project.name, "_no-multiplets_seurat.rds"))

  if(SAC_HTO == ""){
    SAC_HTO = NULL
  }

  # Remove SAC
  Idents(singlets_seurat) <- "HTO_Assignment"
  singlets_seurat <- subset(singlets_seurat, idents = c(SAC_HTO, "Negative"), invert = TRUE)

  # Return all Seurat metadata with barcode
  write.csv(singlets_seurat@meta.data, file = paste0("Processing_QC/", singlets_seurat@project.name, "_metadata_with_bc.csv"), row.names = FALSE)

  # Remove duplicated metadata
  singlets_seurat[["HTO_classification.global"]] <- NULL
  singlets_seurat[["hash.ID"]] <- NULL

  saveRDS(singlets_seurat, paste0("Seurat_objects/", demux_seurat@project.name, "_no-spike-in_singlets_seurat.rds"))

  singlets_seurat
}
