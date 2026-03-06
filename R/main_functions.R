#' Pre-processing Of Raw Data by Project
#'
#' Performs pre-processing of CellRanger data, including background noise removal, HTO demultiplexing and multiplets filtering. Returns a Seurat object.
#' @param seurat_list Seurat list output from merge_pools function.
#' @param isotype Name of the isotype used as ADT control. If no control, FALSE. Default to FALSE.
#' @param denoise_bg Denoising with background mean only. Only used when isotype=FALSE. Default to TRUE.
#' @param project_name (optional) Project name used in the creation of the Seurat object.
#' @param pool_HTO Vector of HTOs in the assay to demultiplex. Default all HTOs.
#' @param prot_bg_threshold Vector of protein size thresholds to define drops as background. Default c(1.5, 3).
#' @param rna_bg_threshold Vector of RNA size thresholds to define drops as background. Default c(1, 2.5).
#' @param method HTO demultiplexing method, HTODemux, MULTI or demuxmix. Default demuxmix.
#' @param make_plots Logical. Default TRUE
#' @param model Character specifying the type of mixture model to used in demuxmix. Either "naive", "reg", "regpos", or "auto" (default).
#' @param clusterInit Optional list of numeric vectors to manually specify the droplet to component assignment used to initialize the EM algorithm of demuxmix. The name of each list element must match a row name of hto. The length of each element must match the number of columns of hto. Only the values 1 and 2 are allowed, where 1 belongs to the negative.
#' @keywords preprocessing
#' @export
#' @import Seurat
#' @import dsb
#' @import Biostrings
#' @import tidyverse
#'

preprocessing <- function(seurat_list, isotype=FALSE, pool_HTO=NULL, method="demuxmix", bg_threshold=seurat_list$bg_thresholds,
                          make_plots=TRUE, denoise_bg=TRUE, project_name=NULL, model="auto", clusterInit=NULL){

  print(paste("Start preprocessing file", project_name))

  # Flag if data includes ADT expression, FALSE if only HTO
  adt.exp <- any(grepl("ADT", names(seurat_list[["filter_data"]]@assays)))
  # Flag if data includes HTO expression, FALSE if not
  hto.exp <- any(grepl("HTO", names(seurat_list[["filter_data"]]@assays)))

  # +/- RNA + ADT + HTO
  if(adt.exp & hto.exp){
    dsb_Seurat <- bg_correction_merged(seurat_list=seurat_list, isotype=isotype, denoise_bg=denoise_bg, pool_HTO=pool_HTO, project_name=project_name,
                                       prot_bg_threshold=bg_threshold[c("prot.low", "prot.high")], rna_bg_threshold=bg_threshold[c("rna.low", "rna.high")])
    demux_Seurat <- demultiplex(dsb_Seurat=dsb_Seurat, pool_HTO=pool_HTO, method=method, mt_threshold=seurat_list[["mt_threshold"]], make_plots=make_plots, model=model, clusterInit=clusterInit)
  }
  # +/- RNA + ADT
  if(adt.exp & !hto.exp){
    demux_Seurat <- bg_correction_merged(seurat_list=seurat_list, isotype=isotype, denoise_bg=denoise_bg, pool_HTO=pool_HTO, project_name=project_name,
                                       prot_bg_threshold=bg_threshold[c("prot.low", "prot.high")], rna_bg_threshold=bg_threshold[c("rna.low", "rna.high")])
  }
  # +/- RNA + HTO
  if(!adt.exp & hto.exp){
    demux_Seurat <- demultiplex(dsb_Seurat=seurat_list[["filter_data"]], pool_HTO=pool_HTO, method=method, mt_threshold=seurat_list[["mt_threshold"]], make_plots=make_plots, model=model, clusterInit=clusterInit)
  }

  write.table(ncol(demux_Seurat), file = paste0("Processing_QC/", project_name, "_rna_qc_cells.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

  print(paste("Finished preprocessing file", project_name))

  demux_Seurat

}


#' Pre-processing Of Raw Data by Pool
#'
#' Performs pre-processing of CellRanger data, including background noise removal, HTO demultiplexing and multiplets filtering. Returns a Seurat object.
#' @param metrics_path Path to CellRanger metrics_summary.
#' @param raw_path Path to CellRanger outs/multi/count/raw_feature_bc_matrix.
#' @param filter_path Path to CellRanger outs/per_sample_outs/SampleName/count/sample_feature_bc_matrix.
#' @param isotype Name of the isotype used as ADT control. Default to HSA. If no control, FALSE.
#' @param denoise_bg Denoising with background mean only. Default to FALSE. Only used when isotype=FALSE.
#' @param project_name (optional) Project name used in the creation of the Seurat object.
#' @param pool_HTO Vector of HTOs in the assay to demultiplex. Default all HTOs.
#' @param mt_threshold Proportion of mitochondrial genes use as cutoff. Default 0.2.
#' @param prot.bg.threshold Vector of protein size thresholds to define drops as background. Default c(1.5, 3).
#' @param rna.bg.threshold Vector of RNA size thresholds to define drops as background. Default c(0, 2.5).
#' @param method Demultiplex method, HTODemux or MULTI. Default HTODemux.
#' @param make_plots Logical. Default TRUE
#' @param to_file Logical. Save metrics to file. Default TRUE
#' @keywords preprocessing
#' @export
#' @import Seurat
#' @import dsb
#' @import Biostrings
#' @import tidyverse
#'

preprocessing_byPool <- function(metrics_path, raw_path, filter_path, batch, isotype="HSA", pool_HTO=NULL, method = "HTODemux", mt_threshold = 0.2, prot.bg.threshold=c(1.5, 3), rna.bg.threshold=c(0, 2.5), make_plots=TRUE, denoise_bg=FALSE, project_name=NULL, to_file=TRUE){

  print(paste("Start preprocessing file", project_name))
  if(to_file){
    extract_metrics(metrics_summary=metrics_path, to_file=to_file, batch = batch)
  }
  
  raw_data <- Read10X(raw_path)
  filter_data <- Read10X(filter_path)

  # Flag if data includes ADT expression, FALSE if only HTO
  adt.exp <- any(!grepl("HTO", rownames(filter_data$`Antibody Capture`)))

  # Check if data has RNA and ADT assays. If only HTO, no correction
  if(length(filter_data) > 1){
    if(adt.exp){ # RNA+ADT+HTO
      dsb_Seurat <- bg_correction(raw_data, filter_data, isotype, denoise_bg, pool_HTO, project_name, mt_threshold, prot.bg.threshold, rna.bg.threshold)
    }else{ # RNA+HTO
      dsb_Seurat <- rna_correction(raw_data, filter_data, pool_HTO, project_name, mt_threshold)
    }
  }else{
    if(adt.exp){ # ADT+HTO
      dsb_Seurat <- adt_correction(raw_data, filter_data, isotype, denoise_bg, pool_HTO, project_name, prot.bg.threshold)
    }else{ # HTO
      dsb_Seurat <- CreateSeuratObject(counts = filter_data, project = project_name, assay = "HTO", min.cells = 0, min.features = 0)
    }
  }

  demux_Seurat <- demultiplex(dsb_Seurat=dsb_Seurat, pool_HTO=pool_HTO, method=method, mt_threshold=mt_threshold, make_plots=make_plots)

  demux_Seurat[["e_barcode"]]<- rownames(demux_Seurat@meta.data)

  # Make scatter plots of ADT in this function?
  print(paste("Finished preprocessing file", project_name))

  demux_Seurat

}

#' Create Input for Change-O
#'
#' Format the merged IMGT VQuest and enclone output for Change-O. Returns a list of data frames (also saved as RData), saves files to tsv/csv and FASTA to use with Change-O.
#' @param summary_path Path to 1_Summary.txt.
#' @param airr_path Path to vquest_airr.tsv.
#' @param enclones List of data frame after top_contig.
#' @param project_name Name of the project.
#' @keywords clone
#' @export
#'

changeo_in <- function(summary_path, airr_path, enclones, project_name){

  v_e <- merge_Vquest(summary_path, airr_path, enclones, project_name)
  v_enclone_to_changeo(v_e, project_name)

  print("Next step, run Change-O with 'filtered_contigs' files.")

  return(v_e)

}

