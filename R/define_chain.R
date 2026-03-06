#' Determine Chains Of enclone Output
#'
#' Determination of chains based on highest read count. It returns a data frame.
#' @param enclone_out enclone output as data.frame, best if it is output from chain_call.
#' @param igm.igd Logical. Call IgM+IgD+ instead of Doublet. Default FALSE.
#' @param project_name Name of the project.
#' @keywords pairing
#' @export
#'

define_chain <- function(enclone_out, igm.igd = FALSE, project_name){

  for(chain in c("HC", "KC", "LC")){

    scoring <- enclone_out %>%
      # Check all possible combinations of chain matching across the 3 seqs
      # When a match is found, divide the read counts for those 2 seqs
      # Invert read counts that are over 1 (this accounts for dividing the higher read count by the lower one)
      mutate(score1 = case_when(chain_1 == chain & chain_2 == chain ~ r_cell1/r_cell2)) %>%
      mutate(score1 = ifelse(score1 >= 1, 1/score1, score1)) %>%

      mutate(score2 = case_when(chain_1 == chain & chain_3 == chain ~ r_cell1/r_cell3)) %>%
      mutate(score2 = ifelse(score2 >= 1, 1/score2, score2)) %>%

      mutate(score3 = case_when(chain_2 == chain & chain_3 == chain ~ r_cell2/r_cell3)) %>%
      mutate(score3 = ifelse(score3 >= 1, 1/score3, score3)) %>%

      # With scores calculated, determine doublets with scores >= 0.5
      # All other scores, determine which read count is higher and put that read count title in the "high_count" column
      mutate(high_count = case_when(score1 < 0.5 ~ colnames(.[,c("r_cell1", "r_cell2")])[max.col(.[,c("r_cell1", "r_cell2")], ties.method = "first")],
                                      score1 >= 0.5 ~ "Doublet",
                                      score2 < 0.5 ~ colnames(.[,c("r_cell1", "r_cell3")])[max.col(.[,c("r_cell1", "r_cell3")], ties.method = "first")],
                                      score2 >= 0.5 ~ "Doublet",
                                      score3 < 0.5 ~ colnames(.[,c("r_cell2", "r_cell3")])[max.col(.[,c("r_cell2", "r_cell3")], ties.method = "first")],
                                      score3 >= 0.5 ~ "Doublet")) %>%

      # Extract the seq# of the read count and combine that with "v_name" to determine which v_name# should be called
      mutate(which_v = case_when(str_detect(high_count, "r_cell") ~ paste("v_name",str_extract(high_count, '\\d'), sep = ""))) %>%

      # Extract the seq# of the read count and combine that with "const" to determine which const# should be called
      mutate(which_c = case_when(str_detect(high_count, "r_cell") ~ paste("const",str_extract(high_count, '\\d'), sep = ""))) %>%

      # Match the v_name# in the which_v column to the column of the same name and use its value as the V gene
      mutate(v = subset(., select = -c(which_v))[cbind(1:nrow(.),match(.$which_v, colnames(.) ))]) %>%

      # Match the const# in the which_c column to the column of the same name and use its value as the C gene
      mutate(c = subset(., select = -c(which_c))[cbind(1:nrow(.),match(.$which_c, colnames(.) ))]) %>%

      # Fill in genes that were doublets
      mutate(v = ifelse(high_count == "Doublet", "Doublet", v)) %>%
      mutate(c = ifelse(high_count == "Doublet", "Doublet", c)) %>%

      # For any instances of a single chain, fill in its V gene to v
      mutate(v = ifelse(chain_1 == chain & c(chain_2 != chain | is.na(chain_2)) & c(chain_3 != chain | is.na(chain_3)), v_name1, v)) %>%
      mutate(v = ifelse(chain_2 == chain & c(chain_1 != chain | is.na(chain_1)) & c(chain_3 != chain | is.na(chain_3)), v_name2, v)) %>%
      mutate(v = ifelse(chain_3 == chain & c(chain_2 != chain | is.na(chain_2)) & c(chain_1 != chain | is.na(chain_1)), v_name3, v)) %>%

      # For any instances of a single chain, fill in its C gene to c
      mutate(c = ifelse(chain_1 == chain & c(chain_2 != chain | is.na(chain_2)) & c(chain_3 != chain | is.na(chain_3)), const1, c)) %>%
      mutate(c = ifelse(chain_2 == chain & c(chain_1 != chain | is.na(chain_1)) & c(chain_3 != chain | is.na(chain_3)), const2, c)) %>%
      mutate(c = ifelse(chain_3 == chain & c(chain_2 != chain | is.na(chain_2)) & c(chain_1 != chain | is.na(chain_1)), const3, c)) %>%
      select(v, c, which_v)

    col_names <- c(paste0("v_", chain), paste0("c_", chain), paste0("which_", chain, "v"))
    colnames(scoring) <- col_names

    enclone_out <- cbind(enclone_out, scoring) # TODO: or join by cell_id?
  }

  # Add this section to proSCessoR::define_chain and add new argument igm.idg = FALSE to pair_analysis
  if(igm.igd){
    # For IgM+IgD+ cells
    # Select cells with IgM and IgD contigs
    tmp <- enclone_out
    tmp$index <- 1:nrow(tmp)
    cols <- grep("const[123]", colnames(tmp), fixed = FALSE)
    cells <- which(rowSums(sapply(tmp[cols], function(x) grepl("IGHD", x))) >= 1)
    enclone_ss <- tmp[cells,]
    cells <- enclone_ss$index[rowSums(sapply(enclone_ss[cols], function(x) grepl("IGHM", x))) >= 1]
    enclone_ss <- tmp[cells,]
    
    # Check if v_name, j_name, cdr3_aa are the same
    res <- enclone_ss %>%
      select(barcode, index, v_name1, v_name2, v_name3, j_name1, j_name2, j_name3, cdr3_aa1, cdr3_aa2, cdr3_aa3) %>%
      rowwise() %>%
      mutate(same.v = ifelse(n_distinct(c_across(c(v_name1, v_name2, v_name3))) <= 2, TRUE, FALSE),
            same.j = ifelse(n_distinct(c_across(c(j_name1, j_name2, j_name3))) <= 2, TRUE, FALSE),
            # TODO: Add check that contig1 and contig2 are always HC
            # Cutoff of 3 was set because of histogram and G001 paper (cutoff for defining clonotypes)
            same.cdr = ifelse(stringdist::stringdist(cdr3_aa1, cdr3_aa2, method = "lv") <= 3, TRUE, FALSE),
            d.pos = ifelse(same.v + same.j + same.cdr == 3, TRUE, FALSE)
            )
    # Save histogram to check if cutoff = 3 still valid
    png(paste0(project_name, "_HCDR3_Lv_dist_IgMIgD.png"))
    hist(stringdist::stringdist(enclone_ss$cdr3_aa1, enclone_ss$cdr3_aa2, method = "lv"), breaks = 20, main = project_name, xlab = "Levenshtein distance between HCDR3")
    dev.off()

    cells <- res$index[res$d.pos]

    enclone_out[cells, "c_HC"] <- "IgM.IgD"
    enclone_out[cells, "v_HC"] <- res[res$d.pos, "v_name1"]
  }

  enclone_out

}
