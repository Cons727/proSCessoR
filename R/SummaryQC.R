#' QC Summary Plots
#'
#' This function generates QC summary plots (all batches) from CellRanger metrics reports and singlets/doublets count distributions.
#' @param batches Vector of batches with batch_ReadMetrics.csv files.
#' @keywords QC
#' @export
#'


SummaryQC <- function(batches){
  # CellRanger metrics report
  # batches <- unique(merged_seurat@meta.data[["orig.ident"]])

  metrics <- list()
  for(i in batches){
    metrics[[i]] <- read.csv(paste0("Processing_QC/", i, "_ReadMetrics.csv"))
    metrics[[i]]$Batch_Pool <- i
  }
  metrics <- do.call(rbind, metrics)
  write.csv(metrics, "Processing_QC/ReadMetrics_Summary.csv")

  lib.type <- unique(metrics$Library.Type)
  metrics$Metric.Name <- gsub("_", " ", metrics$Metric.Name)
  for(i in lib.type){
    ggplot(filter(metrics, Library.Type == i), aes(x = Batch_Pool, y = Metric.Value, color = Batch_Pool)) +
      geom_point(size = 5) +
      facet_wrap(~Metric.Name, scale="free") +
      theme_minimal() +
      theme(axis.text.x=element_blank(), text = element_text(size = 20), legend.position="bottom") +
      labs(x = element_blank(), y = element_blank(), title = i)
    ggsave(paste0("Processing_QC/", i, "_metrics-library_Summary.jpg"), width = 16, height = 6)
  }
}
