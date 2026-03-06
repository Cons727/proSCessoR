#' Get 10x BCR Barcodes
#'
#' Read 10x BCR data and extract the barcodes of droplets classified as cells by Cell Ranger. Files are saved in Processing_QC as txt. Proceed to Pairing Analysis.
#' @param batches Vector of batches to process.
#' @keywords preprocessing
#' @export
#'

BCR_barcodes <- function(batches){
    invisible(ifelse(!dir.exists(file.path(".", "Processing_QC")), dir.create(file.path(".", "Processing_QC")), FALSE))
    
    for(batch in batches){
        path <- paste0(batch, "/per_sample_outs/", batch, "/vdj_b/cell_barcodes.json")
        file <- get_CR_file(file = "cell_barcodes.json", path = path)
        singlets <- rjson::fromJSON(file = file)
        write.table(singlets, file = paste0("Processing_QC/", batch,"_singlets_barcodes.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
    }
}