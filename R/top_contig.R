#' Create FASTA Of Top Contig From enclone Output
#'
#' Select top contig after pairing and generate FASTA file for all chains. It returns a FASTA file and a list of data frames.
#' @param enclone_out enclone output as data.frame, output from pairing.
#' @param project_name Name of project.
#' @keywords pairing
#' @export
#'

top_contig <- function(enclone_out, project_name){

  stringset <- c()
  enclones <- list()

  for(chain in c("HC", "KC", "LC")){

    print(paste0("Analysing chain ", chain))

    which_v <- paste0("which_", chain, "v")
    v <- paste0("v_", chain)

    enclones[[chain]] <- enclone_out %>%
      mutate(Top_contig = case_when(is.na(enclone_out[, which_v]) & (chain_1 != chain | is.na(chain_1)) & (chain_2 != chain | is.na(chain_2)) & (chain_3 != chain | is.na(chain_3)) ~ as.character(0),
                                    is.na(enclone_out[, which_v]) & chain_1 == chain & (chain_2 != chain | is.na(chain_2)) & (chain_3 != chain | is.na(chain_3)) ~ as.character(1),
                                    is.na(enclone_out[, which_v]) & (chain_1 != chain | is.na(chain_1)) & chain_2 == chain & (chain_3 != chain | is.na(chain_3)) ~ as.character(2),
                                    is.na(enclone_out[, which_v]) & (chain_1 != chain | is.na(chain_1)) & (chain_2 != chain | is.na(chain_2)) & chain_3 == chain ~ as.character(3),
                                    str_detect(enclone_out[, which_v], "v_name\\d") == TRUE ~ str_extract(enclone_out[, which_v], '\\d'),
                                    enclone_out[, v] == "Doublet" ~ as.character(0))) %>%
      mutate(full_seq_to_call = paste0("seq", Top_contig)) %>%
      mutate(full_seq = subset(., select = -c(full_seq_to_call))[cbind(1:nrow(.), match(.$full_seq_to_call, colnames(.) ))]) %>%
      mutate(FASTAname = paste0(barcode, "_contig", Top_contig, "_", chain)) %>%

    # vdj_seq_[chain] is not present in final output
      mutate(seq_to_call = paste0("vj_seq", Top_contig)) %>%
      mutate(vdj_seq = subset(., select = -c(seq_to_call))[cbind(1:nrow(.), match(.$seq_to_call, colnames(.)))]) %>%

      mutate(const_to_call = paste("const", Top_contig, sep = "")) %>%
      mutate(c = subset(., select = -c(const_to_call))[cbind(1:nrow(.), match(.$const_to_call, colnames(.) ))]) %>%

    #Top_contig == "1" | Top_contig == "2" | Top_contig == "3" same as Top_HC_contig != "0"?
      mutate(Top_read_count = case_when(Top_contig != "0" ~ paste0("r_cell", Top_contig))) %>%
      mutate(Top_read_count = subset(., select = -c(Top_read_count))[cbind(1:nrow(.), match(.$Top_read_count, colnames(.) ))])%>%

      mutate(Top_UMI_count = case_when(Top_contig != "0" ~ paste0("u_cell", Top_contig))) %>%
      mutate(Top_UMI_count = subset(., select = -c(Top_UMI_count))[cbind(1:nrow(.), match(.$Top_UMI_count, colnames(.) ))])

    fasta_string <- DNAStringSet(enclones[[chain]]$full_seq[!is.na(enclones[[chain]]$full_seq)])
    names(fasta_string) <- enclones[[chain]]$FASTAname[!is.na(enclones[[chain]]$full_seq)]
    stringset <- c(stringset, fasta_string)
  }
  # Change to better concat function above?
  stringset <- do.call(c, stringset)
  writeXStringSet(stringset, paste0("Processing_QC/", project_name, "_TopCalls.fasta"))

  enclones
}
