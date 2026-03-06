#' Extract Metrics From CellRanger Output
#'
#' This function allows you to extract Estimated number of cells, Mean reads per cell, Mean used reads per cell, Antibody reads usable per cell, Mean used reads per cell, Median UMI counts per cell from CellRanger output metrics_summary.csv and can save it as a csv file or print it to console. web_summary.html files will be copied to a Web_Summaries folder in the current directory.
#' @param metrics_summary Path to CellRanger output metrics_summary.csv
#' @param batch Batch name.
#' @param to_file Logical. Save metrics output as a csv file. Default to TRUE.
#' @keywords metrics
#' @export
#'

extract_metrics <- function(metrics_summary, batch, to_file=TRUE){

  # mkdir -p Web_Summaries and copy file to directory
  invisible(ifelse(!dir.exists(file.path(".", "Web_Summaries")), dir.create(file.path(".", "Web_Summaries")), FALSE))
  web_summary <- gsub("metrics_summary.csv", "web_summary.html", metrics_summary, fixed = TRUE)
  file.copy(from = web_summary, to = "Web_Summaries")
  file.rename(from = paste0("Web_Summaries/", basename(web_summary)), to = paste0("Web_Summaries/", batch, "_", basename(web_summary)))

  ReadMetrics <- read.csv(file = metrics_summary, stringsAsFactors = FALSE)

  metrics <- c("Cells",
            "VDJ cells", # CRv10 only
            "Estimated number of cells", # CRv8 only
            "Median UMI counts per cell", "Mean reads per cell",
            "Median reads per cell", # CRv8 only
            "Mean antibody reads usable per cell", "Mean used reads per cell")

  ReadMetrics <- ReadMetrics %>%
    filter(Metric.Name %in% metrics) %>%
    mutate(Metric.Value = as.numeric(gsub(",", "", Metric.Value)),
          Library.Type = case_when(Library.Type == "Antibody Capture" ~ "CITEseq Abs",
                                   Library.Type == "Gene Expression" ~ "GEX",
                                   Library.Type == "VDJ B" ~ "BCR",
                                   Library.Type == "VDJ T" ~ "TCR"),
          Metric.Name = case_when(Metric.Name %in% c("Cells", "VDJ cells", "Estimated number of cells") ~ "number_of_cells",
                                  Metric.Name %in% c("Mean reads per cell", "Median reads per cell") & Library.Type == "GEX" & Category == "Cells" ~ "median_reads_per_cell",
                                  Metric.Name == "Mean reads per cell" & Library.Type == "GEX" ~ "mean_reads_per_cell",
                                  Metric.Name == "Mean used reads per cell" ~ "mean_used_reads_per_cell",
                                  Metric.Name == "Mean antibody reads usable per cell" ~ "mean_usable_reads_per_cell",
                                  Metric.Name == "Median UMI counts per cell" ~ "median_UMI_counts_per_cell")) %>%
    filter(!is.na(Metric.Name)) %>%
    distinct(Library.Type, Metric.Name, Metric.Value, .keep_all = TRUE) %>%
    arrange(desc(Library.Type), Metric.Name) %>%
    mutate(RowNames = paste(Library.Type, Metric.Name, sep = "_")) %>%
    column_to_rownames(var = "RowNames")

  ## Choose to save the data or output it to the console ##
  invisible(ifelse(!dir.exists(file.path(".", "Processing_QC")), dir.create(file.path(".", "Processing_QC")), FALSE))
  if(to_file){
    write.csv(ReadMetrics, file = paste0("Processing_QC/", batch, "_ReadMetrics.csv"), quote = FALSE)
  }
  else {
    print(ReadMetrics)
  }

  # Read depth test
  depth <- ReadMetrics %>%
      filter(grepl("reads", Metric.Name)) %>%
      select(Library.Type, Metric.Value)
  r.depth <- depth$Metric.Value
  names(r.depth) <- depth$Library.Type

  if("GEX" %in% names(r.depth)){
    if(r.depth["GEX"] <= 10000){
      warning(paste0("STOP analysis! \nRead depth not achieved in GEX library! \nPool: ", batch, "\nRead depth = ", r.depth["GEX"]), immediate. = TRUE)
    }
  }
  if("CITEseq Abs" %in% names(r.depth)){
    if(r.depth["CITEseq Abs"] <= 5000){
      warning(paste0("STOP analysis! \nRead depth not achieved in CITEseq Abs library! \nPool: ", batch, "\nRead depth = ", r.depth["CITEseq Abs"]), immediate. = TRUE)
    }
  }
  if("BCR" %in% names(r.depth)){
    if(r.depth["BCR"] <= 10000){
      warning(paste0("STOP analysis! \nRead depth not achieved in BCR library! \nPool: ", batch, "\nRead depth = ", r.depth["BCR"]), immediate. = TRUE)
    }
  }
  if("TCR" %in% names(r.depth)){
    if(r.depth["TCR"] <= 10000){
      warning(paste0("STOP analysis! \nRead depth not achieved in TCR library! \nPool: ", batch, "\nRead depth = ", r.depth["TCR"]), immediate. = TRUE)
    }
  }

}
