process merging_pools {
    container "${params.container__R}"
    // Publish plots and tables
    publishDir "${params.output_directory}", mode: 'copy', overwrite: true

    input:
    // The working directory will contain:
    // P1/* P2/* P3/* .../* (cellranger multi out of each GEM pool)
    path "*"
    path metadata
    path code_repository
    // path to proSCessoR instead of an R script to load the library? Install proSCessoR in the container?

    output:
    // All the processing files are going to be saved in Processing_QC
    path "*"

    script:
    """#!/usr/bin/env Rscript
    print("loading R library")

    devtools::load_all("${code_repository}")
    library(Seurat)
    library(tidyverse)

    print("Merging pools")
    print(list.files(path="."))

    # Set arguments
    project.name <- "${params.project_name}" # Study.ID_batch.ID
    
    SAC.HTO <- c("${params.spikein}")
    if(SAC.HTO == "NULL"){
        SAC.HTO <- NULL
    }

    igm.igd <- FALSE

    batches <- grep(pattern = "P[[:digit:]]", x = list.dirs(full.names = FALSE, recursive = FALSE), value = TRUE)
    print(paste("Processing: ", batches))
    seurat.list <- merge_pools(batches = batches, mt_threshold = ${params.mt_threshold}, project_name = project.name, extract_metrics = TRUE)

    merged.seurat <- preprocessing(seurat_list = seurat.list,
                          method = "demuxmix",
                          isotype = FALSE,
                          project_name = project.name)
    
    invisible(ifelse(!dir.exists(file.path(".", "Seurat_objects")), dir.create(file.path(".", "Seurat_objects")), FALSE))
    saveRDS(merged.seurat, paste0("Seurat_objects/", project.name, "_dsb_whole_dataset_seurat.rds"))

    SummaryQC(batches)

    print("${metadata}")

    singlets.seurat <- join_md(merged.seurat, md_path = "${metadata}", SAC_HTO = SAC.HTO, print_res = TRUE)
    """
}

workflow merge_pools {

    take:
        // Channel containing the output files produced by CellRanger
        input_directory
        metadata
        code_repository

    main:

        // Run merging_pools
        merging_pools(
            input_directory,
            metadata,
            code_repository
        )

    emit:
        merging_pools.out

}
