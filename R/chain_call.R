#' Define Type of Chain and Frequency
#'
#' Defining which type of chain (heavy chain, light kappa chain or light lambda chain) the call was and the frequency each chain was call for each cell. It returns a data frame.
#' @param enclone_out Path to enclone output csv.
#' @keywords pairing
#' @export
#'

chain_call <- function(enclone_out){
  enclone_out <- read.csv(enclone_out, header = TRUE)

  # Type of chain
  enclone_out <- enclone_out %>%
    mutate(chain_1 = case_when(str_detect(v_name1, "IGH") ~ "HC",
                               str_detect(v_name1, "IGK") ~ "KC",
                               str_detect(v_name1, "IGL") ~ "LC",)) %>%
    mutate(chain_2 = case_when(str_detect(v_name2, "IGH") ~ "HC",
                               str_detect(v_name2, "IGK") ~ "KC",
                               str_detect(v_name2, "IGL") ~ "LC",)) %>%
    mutate(chain_3 = case_when(str_detect(v_name3, "IGH") ~ "HC",
                               str_detect(v_name3, "IGK") ~ "KC",
                               str_detect(v_name3, "IGL") ~ "LC",))
  # Frequency of chain
  enclone_out$HC_count <- apply(enclone_out[, c('chain_1', 'chain_2', 'chain_3')], 1, function(x) length(which(x=="HC")))
  enclone_out$KC_count <- apply(enclone_out[, c('chain_1', 'chain_2', 'chain_3')], 1, function(x) length(which(x=="KC")))
  enclone_out$LC_count <- apply(enclone_out[, c('chain_1', 'chain_2', 'chain_3')], 1, function(x) length(which(x=="LC")))

  enclone_out

}
