#' Replaces NA to No or 0
#'
#' Replaces NA to "No" in Single_HC and Light_Chain, and NA to 0 in Paired. It returns a data frame.
#' @param enclone_out Data frame after top_contig.
#' @keywords pairing
#' @export
#'

replace_NA <- function(enclone_out){
  for(chain in c("HC", "KC", "LC")){
    enclone_out[[chain]] <- enclone_out[[chain]] %>%
      mutate(Single_HC = ifelse(is.na(Single_HC), "No", Single_HC)) %>%
      mutate(Light_Chain = ifelse(is.na(Light_Chain), "No", Light_Chain)) %>%
      mutate(Paired = ifelse(is.na(Paired), 0, Paired))
  }

  enclone_out
}
