#' Format V-quest/Enclone output for Change-O
#'
#' Create CellRanger-like FASTA and contig annotation files for use with Change-O. File should not contain barcode overlaps. It returns FASTA file and csv.
#' @param v_enclones List of data frames from merge_Vquest.
#' @param project_name Name of the project.
#' @keywords changeo
#' @export
#'

v_enclone_to_changeo <- function(v_enclones, project_name){

  for(chain in c("HC", "KC", "LC")){

    print(paste0("Analysing chain ", chain))

    v_enclones[[chain]] <- v_enclones[[chain]] %>%
      filter(Paired > 0) %>%
      # Match the filtered_contig_annotation file
      dplyr::rename(v_gene = vGene, d_gene = d_name, j_gene = j_name, c_gene = c, fwr1 = fwr1_aa, fwr1_nt = fwr1, fwr2 = fwr2_aa, fwr2_nt = fwr2, fwr3 = fwr3_aa, fwr3_nt = fwr3, fwr4 = fwr4_aa, fwr4_nt = fwr4, cdr1 = cdr1_aa, cdr1_nt = cdr1, cdr2 = cdr2_aa, cdr2_nt = cdr2, cdr3 = cdr3_aa, cdr3_nt = cdr3, contig_sequence = sequence, productive = productive, full_length = complete_vdj, reads = Top_read_count, umis = Top_UMI_count) %>%
      mutate(length = nchar(contig_sequence)) %>%
      mutate(chain = paste0("IG", strsplit(chain, split=NULL)[[1]][1])) %>%
      mutate(is_cell = "TRUE") %>%
      mutate(high_confidence = "TRUE") %>%
      mutate(contig_id = paste(barcode, "contig", Top_contig, sep = "_")) %>%
      select("barcode","is_cell","contig_id","high_confidence","length","chain","v_gene","d_gene","j_gene","c_gene","full_length","productive","fwr1","fwr1_nt","cdr1","cdr1_nt","fwr2","fwr2_nt","cdr2","cdr2_nt","fwr3","fwr3_nt","cdr3","cdr3_nt","fwr4","fwr4_nt","reads","umis","contig_sequence") %>%
      # Remove any rows where there was no v_gene call
      filter(!is.na(v_gene))
  }

  topcalls <- do.call(bind_rows, v_enclones) %>%
    arrange(barcode)

  topcalls[is.na(topcalls)] <- ""

  topcalls_dna <- DNAStringSet(topcalls$contig_sequence)
  names(topcalls_dna) <- topcalls$contig_id

  topcalls %>% select(-contig_sequence) %>%
    write.csv(file = paste0("Processing_QC/", project_name, "_post_e-Vq_filtered_contig_annotations.csv"), row.names = FALSE, quote = FALSE)

  writeXStringSet(topcalls_dna, paste0("Processing_QC/", project_name, "_post_e-Vq_filtered_contig.fasta"))
}
