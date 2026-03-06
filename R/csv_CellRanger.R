#' Generate Input CSV File For Cell Ranger
#'
#' For each unique "sample name" in the specified directory, a new sample_name_CellRanger_Submission.csv file is created in the current directory.
#' @param input_path Path to the directory where the FastQ files are located; all files will be processed.
#' @param feature_ref_path Path to the feature reference file.
#' @param expect_cells Cells expected in batch.
#' @keywords cellranger
#' @export
#'

# library(magrittr)
# library(dplyr)

csv_CellRanger <- function(input_path, feature_ref_path, expect_cells, gex.ref, vdj.ref){

  files <- list.files(input_path)
  batches <- sapply(strsplit(files, split = "_(GEX|FBC|BCR)_.*", fixed = FALSE), '[', 1)
  seq_name <- sapply(strsplit(files, split = "_S[0-9]*_.*", fixed = FALSE), '[', 1)
  library <- ifelse(grepl("GEX", seq_name), "gene expression",
                    ifelse(grepl("FBC", seq_name), "antibody capture", "vdj-b"))
  libraries <- data.frame(Batch = batches, seqs = seq_name, paths = input_path, lanes = "any", features = library)
  libraries <- libraries[!duplicated(libraries), ]

  for(batch in unique(batches)){
    library <- libraries %>%
      filter(Batch == batch) %>%
      select(-Batch)
    table <- data.frame(seqs = c("[gene-expression]",
                                 "reference",
                                 "expect-cells",
                                 "[feature]",
                                 "reference",
                                 "[vdj]",
                                 "reference",
                                 "[libraries]",
                                 "fastq_id"),
                        paths = c("",
                                 gex.ref,
                                 expect_cells,
                                 "",
                                 feature_ref_path,
                                 "",
                                 vdj.ref,
                                 "",
                                 "fastqs"),
                        lanes = c(rep("", 8), "lanes"),
                        features = c(rep("", 8), "feature_types"))
    table <- rbind(table, library)
    write.table(
      table,
      file = paste0(batch, "_CellRanger_Submission.csv"),
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE,
      sep=","
    )
  }
}

# # Get the command line arguments from the user
# options <- commandArgs(trailingOnly = TRUE)
#
# # Run the script
# csv_CellRanger(
#   options[1],
#   options[2],
#   options[3],
#   options[4],
#   options[5]
# )
