#' Get Cell Ranger File Path
#'
#' Get path to a Cell Ranger file which can be saved in either outs or per_sample_outs directory. This command allows to use CellRanger Cirro output with proSCessoR functions.
#' @param file String. Pattern or filename to retrieve path.
#' @param path String. Path to parent directory of file.
#' @keywords preprocessing
#' @export
#'

get_CR_file <- function(path, file){
    parent.dir <- str_split_i(path, pattern = "/", 1)
    outs.path <- sub(pattern = parent.dir, replacement = paste0(parent.dir, "/outs"), x = path)
    v10.path <- gsub(pattern = "multi|count", replacement = "", x = path)
    v10.outs.path <- gsub(pattern = "multi|count", replacement = "", x = outs.path)

    file <- list.files(path = c(path, outs.path, v10.path, v10.outs.path),
                       pattern = file,
                       full.names = TRUE,
                       recursive = FALSE)
    # print(file)
    return(file)
}