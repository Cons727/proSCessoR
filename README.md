# proSCessoR

R package for preprocessing CITE-seq and V(D)J data

![alt text](https://github.com/Cons727/proSCessoR/blob/main/flowchart-General.v0.2.0.png)

# Installation

Install imported packages

    packages <- c('Matrix', 'airr', 'alakazam', 'data.table', 'dplyr', 'dsb', 'ggplot2', 'gplots', 'lubridate', 'methods', 'purrr', 'readr', 'scoper', 'shazam', 'tibble', 'tidyr')
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages)

proSCessoR pipeline requires [Change-O](https://changeo.readthedocs.io/en/stable/overview.html) for the clonotype analysis.

    pip3 install changeo --user

proSCessoR pipeline requires Bioconductor packages for sequence manipulation.

    if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

    BiocManager::install(c("Biostrings", "GenomicAlignments", "limma", "demuxmix"))

Clone proSCessoR repo

    git clone https://github.com/Cons727/proSCessoR.git

If you want to use the pipeline in Fred Hutch's HPC, make files executable

    chmod +x path/to/proSCessoR/bin/*

The easiest way to install proSCessoR is using devtools. Download the package, make sure to have devtools installed and run

    devtools::install("path/to/proSCessoR")
    # devtools::load_all("path/to/proSCessoR")
    library(proSCessoR)

# Workflow

## FastQC reports
Before running the workflow, make sure to generate the QC reports of the raw FASTQ files from the computing cluster. This step is only done once per every sequencing run (sequencing cell ID). Go to the desired directory and run the bash script from the command line.

    cd output/dir
    sbatch -c 16 --mail-user=cmarinim@fredhutch.org path/to/proSCessoR/sbatch_concat_fastqc.sh /shared/ngs/illumina/230412_A00613_0534_BHJ7J/Unaligned/

The FastQC reports for each sample and one MultiQC report of all files in the directory are saved in FastQC_Reports.

## Processing of 10x data
In R, set the working directory where you want the output to be saved. All the intermidiate results are saved in Processing_QC directory. Set other arguments used throughout the pipeline such as the project name, batch name, the hashtag label of a spike-in control or samples that are not going to be analyzed, and if IgM+IgD+ calling is requiered.

    wkdir <- "SingleCell/HVTN302/"
    setwd(wkdir)

    # Set arguments
    project.name <- "HVTN302_BCS2036"
    batch.name <- "BCS2036"
    SAC.HTO <- NULL
    igm.igd <- FALSE

The first step is merging all the files together. Make sure 'batches' only contains paths to CellRanger output. 'merge_pools' will return a list of Seurat objects for further use, and it will save a plot in the current directory which is employed to select the background area for protein normalization. If 'extract_metrics'=TRUE, then Cell Ranger web summaries will be copied to Web_Summaries.

    batches <- grep("^P[0-9]", list.dirs(path = wkdir, full.names = FALSE, recursive = FALSE), value = TRUE)
    seurat.list <- merge_pools(batches=batches, mt_threshold = 0.2, project_name = "HVTN302", extract_metrics = TRUE)

After merging the GEM pools, the 'preprocessing' function can be run, which includes background correction and HTO demultiplexing. Demultiplexing can be performed using demuxmix ("demuxmix"), HTODemux ("HTODemux"), or MULTIseqDemux ("MULTI").
    
    merged.seurat <- preprocessing(seurat_list = seurat.list,
                          method = "demuxmix",
                          isotype = FALSE,
                          project_name = project.name)
    
    invisible(ifelse(!dir.exists(file.path(".", "Seurat_objects")), dir.create(file.path(".", "Seurat_objects")), FALSE))
    saveRDS(so, paste0("Seurat_objects/", project.name, "_dsb_whole_dataset_seurat.rds"))

If 'extract_metrics'=TRUE was used in preprocessing function, generate quality control plots and tables.

    SummaryQC(batches)

Add the metadata to the merged Seurat object. The metadata file must contain HTO and Batch_Pool information. The following function returns a clean Seurat object (without empty droplets and no spike-in control (SAC.HTO argument)) and saves the QC'ed singlet barcodes per batch (or GEM pool) for pairing and clonotyping steps. If the metadata includes SortPop_Small, 'print_res'=TRUE will return demultiplexing heatmaps with short sample IDs.

    singlets.seurat <- join_md(merged.seurat, md_path = "HVTN302_metadata.csv", SAC_HTO = SAC.HTO, print_res = TRUE)

Make sure the QC plots and tables look good before continuing to the next step.

### Optional: Preprocessing by pool
The first step is preprocessing all the files. Make sure 'batches' only contains paths to CellRanger output. Change the working directory for every batch to avoid overwriting (TODO: add file_name argument to extract_metrics). The 'preprocessing_byPool' function includes the background correction and demultiplexing. It will return a Seurat object with only singlets, and the metrics reports are saved as *ReportMetrics*. After running each batch, merge or integrate all the datasets together.

    batches <- list.dirs(path = wkdir, full.names = FALSE, recursive = FALSE)
    obj_seurat <- list()

    for(batch in batches){
      setwd(paste0(wkdir, batch))
      
      metrics_path <- paste0(wkdir, batch, "/outs/per_sample_outs/", batch, "/metrics_summary.csv")
      raw_path <- paste0(wkdir, batch, "/outs/multi/count/raw_feature_bc_matrix")
      filter_path <- paste0(wkdir, batch, "/outs/per_sample_outs/", batch, "/count/sample_filtered_feature_bc_matrix")
        obj_seurat[[batch]] <- preprocessing_byPool(metrics_path = metrics_path,
                                                raw_path = raw_path,
                                                filter_path = filter_path,
                                                project_name = batch,
                                                method = "demuxmix",
                                                isotype = FALSE)
    }

    saveRDS(obj_seurat, "dsb_full_list_seurats.rds")

    merged.seurat <- merge(obj_seurat[[1]], y = obj_seurat[-1], add.cell.ids = names(obj_seurat))

    saveRDS(merged.seurat, "dsb_full_seurat.rds")

## Pairing and clonotype analyses

The pairing and clonotyping steps are performed in one simple function, which will be run in a subshell of the HPC. Since each python and R versions use different C++ libraries, it is important to specify the python site-packages path and the python module.

    pair_and_clone(outname = project.name,
               ref = "/home/OGRDB_ref/ref-2026.01.26/airr_c_human",
               db_name = "ogrdb",
               vdj = "vdj_b",
               igm_igd = "FALSE",
               python_site_packages = "/home/.local/lib/python3.8/site-packages",
               python_module = "Python/3.8.2-GCCcore-9.3.0")

When the process has finished, the final result is saved as *today_project.name_ChangeO-IgSADIE.tsv* in the working directory.

If you only need the pairing analysis without performing the clonotyping, use the following function

    pair_and_merge(outname = project.name,
               ref = "/home/OGRDB_ref/ref-2026.01.26/airr_c_human",
               db_name = "ogrdb",
               vdj = "vdj_b",
               igm_igd = "FALSE",
               python_site_packages = "/home/.local/lib/python3.8/site-packages")

When the process has finished, the final result is saved as *today_project.name_IgSADIE.tsv* in the working directory.