#' Run enclone
#'
#' Run enclone from files listed in vector batches. The function calls enclone from the command line. Saves merged_enclone_output.csv to current directory.
#' @param batches Directory names of each batch to run.
#' @param classifications_path Path to dsb_sample_classifications.tsv. Default "./project_name_dsb_sample_classifications.tsv".
#' @param project_name Name of the project.
#' @keywords enclone
#' @export
#'

run_enclone <- function(batches, project_name, classifications_path = NULL){

  if(is.null(classifications_path)){
    classifications_path <- paste0("Processing_QC/", project_name, "_dsb_sample_classifications.tsv")
  }

  columns <- " PCOLS=barcode,HTO_Assignment,Batch_Pool,clonotype_id,exact_subclonotype_id,cred,gex,n,datasets,origins,donors,datasets_cell,origins_cell,donors_cell,n_gex,near,far,dref,dref_aa,clonotype_ncells,seq1,v_name1,d_name1,j_name1,const1,d_donor1,vj_seq1,vj_seq_nl1,vj_aa1,vj_aa_nl1,fwr1_dna1,fwr1_aa1,cdr1_dna1,cdr1_aa1,fwr2_dna1,fwr2_aa1,cdr2_dna1,cdr2_aa1,fwr3_dna1,fwr3_aa1,cdr3_dna1,cdr3_aa1,fwr4_dna1,fwr4_aa1,fwr1_len1,fwr2_len1,fwr3_len1,fwr4_len1,d1_name1,d2_name1,d1_score1,d2_score1,d_delta1,d_univ1,d_frame1,d_start1,cdr1_dna_ref1,cdr2_dna_ref1,cdr1_aa_north1,cdr1_aa_ref1,cdr2_aa_north1,cdr2_aa_ref1,cdr3_aa_north1,cdr3_aa_conx1,cdr3_aa_conp1,fwr1_dna_ref1,fwr2_dna_ref1,fwr3_dna_ref1,fwr4_dna_ref1,fwr1_aa_ref1,fwr2_aa_ref1,fwr3_aa_ref1,fwr4_aa_ref1,ulen1,vjlen1,clen1,notes1,u_cell1,r_cell1,var1,comp1,edit1,cigar1,v_id1,d_id1,j_id1,const_id1,utr_id1,utr_name1,cdr3_start1,v_start1,var_indices_dna1,var_indices_aa1,share_indices_dna1,share_indices_aa1,seq2,v_name2,d_name2,j_name2,const2,d_donor2,vj_seq2,vj_seq_nl2,vj_aa2,vj_aa_nl2,fwr1_dna2,fwr1_aa2,cdr1_dna2,cdr1_aa2,fwr2_dna2,fwr2_aa2,cdr2_dna2,cdr2_aa2,fwr3_dna2,fwr3_aa2,cdr3_dna2,cdr3_aa2,fwr4_dna2,fwr4_aa2,fwr1_len2,fwr2_len2,fwr3_len2,fwr4_len2,d1_name2,d2_name2,d1_score2,d2_score2,d_delta2,d_univ2,d_frame2,d_start2,cdr1_dna_ref2,cdr2_dna_ref2,cdr1_aa_north2,cdr1_aa_ref2,cdr2_aa_north2,cdr2_aa_ref2,cdr3_aa_north2,cdr3_aa_conx2,cdr3_aa_conp2,fwr1_dna_ref2,fwr2_dna_ref2,fwr3_dna_ref2,fwr4_dna_ref2,fwr1_aa_ref2,fwr2_aa_ref2,fwr3_aa_ref2,fwr4_aa_ref2,ulen2,vjlen2,clen2,notes2,u_cell2,r_cell2,var2,comp2,edit2,cigar2,v_id2,d_id2,j_id2,const_id2,utr_id2,utr_name2,cdr3_start2,v_start2,var_indices_dna2,var_indices_aa2,share_indices_dna2,share_indices_aa2,seq3,v_name3,d_name3,j_name3,const3,d_donor3,vj_seq3,vj_seq_nl3,vj_aa3,vj_aa_nl3,fwr1_dna3,fwr1_aa3,cdr1_dna3,cdr1_aa3,fwr2_dna3,fwr2_aa3,cdr2_dna3,cdr2_aa3,fwr3_dna3,fwr3_aa3,cdr3_dna3,cdr3_aa3,fwr4_dna3,fwr4_aa3,fwr1_len3,fwr2_len3,fwr3_len3,fwr4_len3,d1_name3,d2_name3,d1_score3,d2_score3,d_delta3,d_univ3,d_frame3,d_start3,cdr1_dna_ref3,cdr2_dna_ref3,cdr1_aa_north3,cdr1_aa_ref3,cdr2_aa_north3,cdr2_aa_ref3,cdr3_aa_north3,cdr3_aa_conx3,cdr3_aa_conp3,fwr1_dna_ref3,fwr2_dna_ref3,fwr3_dna_ref3,fwr4_dna_ref3,fwr1_aa_ref3,fwr2_aa_ref3,fwr3_aa_ref3,fwr4_aa_ref3,ulen3,vjlen3,clen3,notes3,u_cell3,r_cell3,var3,comp3,edit3,cigar3,v_id3,d_id3,j_id3,const_id3,utr_id3,utr_name3,cdr3_start3,v_start3,var_indices_dna3,var_indices_aa3,share_indices_dna3,share_indices_aa3 PCELL"

  for(batch in batches){
    print(paste0("Running enclone for batch ", batch))
    flags <- paste0("BCR_GEX=", batch, "/outs",
                  " BC=", classifications_path,
                  " KEEP_CELL_IF=\"Status == 'Singlet' && Batch_Pool == '", batch, "'\"",
                  " POUT=", batch, "/", batch, "_enclone_output.csv",
                  columns)
    system2(Sys.getenv("ENCLONE"), args = flags, stdout = paste0(batch, "/", batch, "_enclone.out"), stderr = paste0(batch, "/", batch, "_enclone.log"))
  }
  # Concatenate enclone_outs
  enclones <- paste0(batches, "/", batches, "_enclone_output.csv")
  enclones <- do.call(rbind, lapply(enclones, read.csv))

  enclones <- enclones %>%
    dplyr::rename(e_barcode = barcode) %>%
    mutate(barcode = paste(Batch_Pool, e_barcode, sep = "_")) %>%
    relocate(barcode, .after = e_barcode)
  write.csv(enclones, paste0("Processing_QC/", project_name, "_merged_enclone_output.csv"), row.names = FALSE)

}
