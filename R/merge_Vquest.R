#' Merge IMGT HighV-Quest Results with Enclone
#'
#' Merge 1_Summary.txt and vquest_airr.tsv files with enclone data frame. Returns a list of data frames and saves output as tsv and RData.
#' @param summary_path Path to 1_Summary.txt.
#' @param airr_path Path to vquest_airr.tsv.
#' @param enclones List of data frame after top_contig.
#' @param project_name Name of project.
#' @keywords vquest
#' @export
#'

merge_Vquest <- function(summary_path, airr_path, enclones, project_name){

  vquest_summary <- read.table(summary_path, sep = "\t", header = TRUE, fill = TRUE, check.names = FALSE)
  vquest_airr <- read.table(airr_path, sep = "\t", header = TRUE, fill = TRUE, check.names = FALSE)

  #Replace sequence_id and Sequence ID with FastaName from enclones. 50 character limit
  vquest_summary$`Sequence ID` <- c(enclones[[1]]$FASTAname[!is.na(enclones[[1]]$full_seq)],
                                    enclones[[2]]$FASTAname[!is.na(enclones[[2]]$full_seq)],
                                    enclones[[3]]$FASTAname[!is.na(enclones[[3]]$full_seq)])
  vquest_airr$sequence_id <- vquest_summary$`Sequence ID`

  vquest_airr <- vquest_airr %>%
    left_join(select(vquest_summary, `Sequence ID`, `AA JUNCTION`, `CDR1-IMGT length`, `CDR2-IMGT length`, `CDR3-IMGT length`, `CDR-IMGT lengths`,
                     `FR-IMGT lengths`, `JUNCTION frame`, `V-REGION identity nt`, `V-REGION identity % (with ins/del events)`,
                     `V-REGION identity nt (with ins/del events)`, `J-REGION identity nt`, `D-REGION reading frame`, `V-DOMAIN Functionality comment`,
                     `V-REGION potential ins/del`, `J-GENE and allele comment`, `Analysed sequence length`, `Sequence analysis category`),
              by = c("sequence_id" = "Sequence ID"))

  # Split the sequence ID into Barcode and Contig, removing the word "contig" but keeping the original column.
  vquest_airr <- vquest_airr %>%
    mutate(contig = grep("[HKL]C", unlist(strsplit(sequence_id, "_")), value = TRUE)) %>%
    mutate(barcode = gsub("_contig[123]_[HKL]C", "", sequence_id))

  vquest_airr <- vquest_airr %>%
    mutate(chain = case_when(str_detect(v_call, "IGH") ~ "HC",
                             str_detect(v_call, "IGK") ~ "KC",
                             str_detect(v_call, "IGL") ~ "LC")) %>%
    relocate(chain, .after = contig) %>%
    dplyr::rename(vGene_and_allele = v_call) %>%
    mutate(vGene = str_extract(vGene_and_allele, "IG[HKL]V\\d+\\w*\\-\\w*\\d+\\w*\\-*\\w*\\d*")) %>%
    mutate(d_name = str_extract(d_call, "IG[HKL]D\\d+\\-\\d+")) %>%
    mutate(j_name = str_extract(j_call, "IG[HKL]J\\d+"))

  for(chain in c("HC", "KC", "LC")){

    print(paste0("Analysing chain ", chain))

    enclones[[chain]] <- left_join(enclones[[chain]], vquest_airr, by = c("FASTAname" = "sequence_id")) %>%
      select(-c(c_call, contig, chain)) %>%
      dplyr::rename(barcode = barcode.x)

    if(chain != "HC"){
      enclones[[chain]] <- enclones[[chain]] %>%
        select(-c(d_call, d_score, d_identity, d_support, d_cigar, d_sequence_start, d_sequence_end, d_germline_start, d_germline_end, d_alignment_start, d_alignment_end, d_sequence_alignment, d_sequence_alignment_aa, d_germline_alignment, d_germline_alignment_aa, Unpaired, Single_HC, Light_Chain, cred))
    }

    enclones[[chain]] <- enclones[[chain]] %>%
      mutate(Top_pct_v_mutation = 100 - v_identity) %>%
      mutate(Top_CDR3_length = junction_aa_length - 2) %>%
      dplyr::rename(v_identity_pct = v_identity, j_identity_pct = j_identity)

    column_order <- c("barcode","Paired","Unpaired","Single_HC","Light_Chain","cred","vGene","vGene_and_allele","Top_pct_v_mutation","c","Top_CDR3_length","Top_contig","Top_read_count","Top_UMI_count","v_identity_pct","d_name","j_name","j_identity_pct","AA JUNCTION","CDR1-IMGT length","CDR2-IMGT length","CDR3-IMGT length","CDR-IMGT lengths","FR-IMGT lengths","d_call","j_call","sequence","productive","rev_comp","complete_vdj","vj_in_frame","stop_codon","sequence_alignment","sequence_alignment_aa","germline_alignment","germline_alignment_aa","junction","junction_aa","JUNCTION frame","v_score","V-REGION identity nt","V-REGION identity % (with ins/del events)","V-REGION identity nt (with ins/del events)","j_score","J-REGION identity nt","D-REGION reading frame","V-DOMAIN Functionality comment","V-REGION potential ins/del","J-GENE and allele comment","insertions","deletions","5prime_trimmed_n_nb","3prime_trimmed_n_nb","Analysed sequence length","Sequence analysis category","np1","np1_aa","np2","np2_aa","cdr1","cdr1_aa","cdr2","cdr2_aa","cdr3","cdr3_aa","fwr1","fwr1_aa","fwr2","fwr2_aa","fwr3","fwr3_aa","fwr4","fwr4_aa","v_support","v_cigar","d_score","d_identity","d_support","d_cigar","j_support","j_cigar","c_score","c_identity","c_support","c_cigar","v_sequence_start","v_sequence_end","v_germline_start","v_germline_end","v_alignment_start","v_alignment_end","d_sequence_start","d_sequence_end","d_germline_start","d_germline_end","d_alignment_start","d_alignment_end","j_sequence_start","j_sequence_end","j_germline_start","j_germline_end","j_alignment_start","j_alignment_end","cdr1_start","cdr1_end","cdr2_start","cdr2_end","cdr3_start","cdr3_end","fwr1_start","fwr1_end","fwr2_start","fwr2_end","fwr3_start","fwr3_end","fwr4_start","fwr4_end","v_sequence_alignment","v_sequence_alignment_aa","d_sequence_alignment","d_sequence_alignment_aa","j_sequence_alignment","j_sequence_alignment_aa","c_sequence_alignment","c_sequence_alignment_aa","v_germline_alignment","v_germline_alignment_aa","d_germline_alignment","d_germline_alignment_aa","j_germline_alignment","j_germline_alignment_aa","c_germline_alignment","c_germline_alignment_aa","junction_length","junction_aa_length","np1_length","np2_length","n1_length","n2_length","p3v_length","p5d_length","p3d_length","p5j_length","rearrangement_id","repertoire_id","rearrangement_set_id","sequence_analysis_category","d_number","junction_decryption")

    enclones[[chain]] <- enclones[[chain]] %>%
      select(any_of(column_order))
  }

  # Save list without changing names
  v_enclones <- enclones

  for(chain in c("HC", "KC", "LC")){
    colnames(enclones[[chain]])[which(colnames(enclones[[chain]])=="vGene"):ncol(enclones[[chain]])] <- paste(colnames(enclones[[chain]])[which(colnames(enclones[[chain]])=="vGene"):ncol(enclones[[chain]])], chain, sep = "_")
  }

  enclones <- purrr::reduce(enclones, full_join, by = "barcode") %>%
    select(-c(Paired.y, Paired.x))
  colnames(enclones) <- gsub(".x", "", colnames(enclones))
  write.table(enclones, paste0("Processing_QC/", project_name, "_v-enclone.tsv"), quote = F, row.names = F, sep = "\t")

  save(v_enclones, file = paste0("Processing_QC/", project_name, "_v_enclones.RData"))

  return(v_enclones)
}
