#' B Cell Pairing Analysis Of IgBLAST Output
#'
#' Pairing analysis after chain determination. It returns a data frame.
#' @param anno Output from run_IgSADIE as data.frame.
#' @param project_name Project name.
#' @param singlets_prod Path to singlets_prod_contig.tsv.
#' @param kl.light Logical. Keep cells with Kappa and Lambda chains. Default TRUE.
#' @keywords pairing
#' @export
#'

pairing <- function(anno, project_name, singlets_prod = NULL, kl.light = TRUE){
  
  if(is.null(singlets_prod)){
    singlets_prod <- paste0("Processing_QC/", project_name, "_singlets_prod_contig.tsv")
  }
  
  singlets_prod <- read.delim(singlets_prod, sep = "\t")

  # Join anno with singlets_prod_f[cell_id, const, reads, umis, selected_contig] by sequence_id
  # Doublets will have NA in anno data
  anno.pair <- singlets_prod %>%
    select(cell_id, sequence_id, const, reads, umis, selected_contig, locus) %>%
    full_join(anno, ., by = "sequence_id") %>%
    mutate(locus = ifelse(is.na(locus.y), locus.x, locus.y), # keep locus information of all cells (singlets + doublets)
           v_mutation_pct = 100 - v_identity,
           cdr3_aa_length = nchar(cdr3_aa), # same as junction_aa_length - 2
           productive = ifelse(productive == "", "F", productive)) %>%
    select(-locus.x, -locus.y)

  # Test for computer artifacts (different locus in VDJ-C genes)
  # IGH requires 4 matches (V, D, J, C)
  # IGK and IHL require 3 matches (V, J, C)
  anno.pair <- anno.pair %>%
    rowwise() %>%
    mutate(check = ifelse(locus == "IGH", sum(str_detect(string = c(v_call, d_call, j_call, c_call), pattern = paste0("^$|", locus))),
                           sum(str_detect(string = c(v_call, j_call, c_call), pattern = paste0("^$|", locus)))),
           artifact = case_when(locus == "IGH" ~ ifelse(check == 4, FALSE, TRUE),
                               locus %in% c("IGK", "IGL") ~ ifelse(check == 3, FALSE, TRUE)),
          real_prod = ifelse(!liable & !artifact & productive == "T", TRUE, FALSE)) %>%
    select(-check, -v_penalty, -d_penalty, -j_penalty, -allow_vdj_overlap,
          -file, -penalty_rank, -penalties)

  anno.pair.wide <- anno.pair %>%
    pivot_wider(id_cols = cell_id,
                names_from = locus,
                names_glue = "{.value}_{locus}",
                values_from = grep("cell_id", colnames(anno.pair), value = TRUE, invert = TRUE)) %>%
    mutate(Single_HC = case_when(selected_contig_IGH == "Doublet-IGH" ~ "Doublet",
                                !is.na(selected_contig_IGH) & real_prod_IGH ~ "Yes",
                                .default = "No"),
          Light_Chain = case_when(!(is.na(selected_contig_IGK) | grepl("Doublet", selected_contig_IGK)) &
                                    !(is.na(selected_contig_IGL) | grepl("Doublet", selected_contig_IGL)) &
                                    real_prod_IGK & real_prod_IGL ~ "KL",
                                  selected_contig_IGK == "Doublet-IGK" ~ "Doublet-IGK",
                                  selected_contig_IGL == "Doublet-IGL" ~ "Doublet-IGL",
                                  !is.na(selected_contig_IGK) &
                                    locus_IGK == "IGK" & 
                                    real_prod_IGK ~ "K",
                                  !is.na(selected_contig_IGL) &
                                    locus_IGL == "IGL" &
                                    real_prod_IGL ~ "L",
                                  .default = "No"),
          Paired = case_when(Single_HC == "Yes" & Light_Chain == "KL" ~ 2,
                             Single_HC == "Yes" &
                                Light_Chain %in% c("K", "L") ~ 1,
                              .default = 0)
    ) %>%
    relocate(Paired, Single_HC, Light_Chain, const_IGH, const_IGK, const_IGL, .after = cell_id) %>%
    relocate(v_mutation_pct_IGH, v_mutation_pct_IGK, v_mutation_pct_IGL,
            cdr3_aa_length_IGH, cdr3_aa_length_IGK, cdr3_aa_length_IGL, .after = c_call_IGL) %>%
    # If this information is of interest, comment line
    select(-selected_contig_IGH, -selected_contig_IGK, -selected_contig_IGL,
          -locus_IGH, -locus_IGK, -locus_IGL,
          -sequence_id_IGH, -sequence_id_IGK, -sequence_id_IGL,
          -real_prod_IGH, -real_prod_IGK, -real_prod_IGL)
  write.table(anno.pair.wide, paste0("Processing_QC/", project_name,"_sadie_wide.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  if(kl.light){
    paired <- anno.pair.wide$cell_id[anno.pair.wide$Paired > 0]
  } else{
    paired <- anno.pair.wide$cell_id[anno.pair.wide$Paired == 1]
  }

  # Filter unpaired (non-productive, liable and artifacts) from airr for changeo-10x.sh (without AssignGenes.py)(modified 2023.12.08)
  anno.pair %>%
    filter(cell_id %in% paired, real_prod) %>%
    select(-selected_contig, -real_prod) %>%
    write.table(paste0("Processing_QC/", project_name,"_ig.sadie_paired_airr.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

}



#' T Cell Pairing Analysis Of IgBLAST Output
#'
#' Pairing analysis after chain determination. Works with both A/B and G/D. It returns a data frame.
#' @param anno Output from run_IgSADIE as data.frame.
#' @param project_name Project name.
#' @param singlets_prod Path to singlets_prod_contig.tsv.
#' @keywords pairing
#' @export
#'

tcr.pairing <- function(anno, project_name, singlets_prod = NULL){
  
  if(is.null(singlets_prod)){
    singlets_prod <- paste0("Processing_QC/", project_name, "_singlets_prod_contig.tsv")
  }
  
  singlets_prod <- read.delim(singlets_prod, sep = "\t")

  # Join anno with singlets_prod_f[cell_id, const, reads, umis, selected_contig] by sequence_id
  # Doublets will have NA in anno data
  anno.pair <- singlets_prod %>%
    select(cell_id, sequence_id, const, reads, umis, selected_contig, locus) %>%
    full_join(anno, ., by = "sequence_id") %>%
    mutate(locus = ifelse(is.na(locus.y), locus.x, locus.y), # keep locus information of all cells (singlets + doublets)
           v_mutation_pct = 100 - v_identity,
           cdr3_aa_length = nchar(cdr3_aa), # same as junction_aa_length - 2
           productive = ifelse(productive == "", "F", productive)) %>%
    select(-locus.x, -locus.y)

  # Test for computer artifacts (different locus in VDJ-C genes)
  # TRB and TRD require 4 matches (V, D, J, C)
  # TRA and TRG require 3 matches (V, J, C)
  anno.pair <- anno.pair %>%
    rowwise() %>%
    mutate(check = case_when(
             locus %in% c("TRB", "TRD") ~ sum(str_detect(string = c(v_call, d_call, j_call, c_call), pattern = paste0("^$|", locus))),
             locus %in% c("TRA", "TRG") ~ sum(str_detect(string = c(v_call, j_call, c_call), pattern = paste0("^$|", locus)))),
           artifact = case_when(locus %in% c("TRB", "TRD") ~ ifelse(check == 4, FALSE, TRUE),
                                locus %in% c("TRA", "TRG") ~ ifelse(check == 3, FALSE, TRUE)),
          real_prod = ifelse(!liable & !artifact & productive == "T", TRUE, FALSE)) %>%
          ungroup() %>%
    select(-check, -v_penalty, -d_penalty, -j_penalty, -allow_vdj_overlap)

  anno.pair.wide <- anno.pair %>%
    pivot_wider(id_cols = cell_id,
                names_from = locus,
                names_glue = "{.value}_{locus}",
                values_from = grep("cell_id", colnames(anno.pair), value = TRUE, invert = TRUE)) %>%
    mutate(Single_HC = case_when(selected_contig_TRB == "Doublet-TRB" ~ "Doublet-TRB",
                                  selected_contig_TRD == "Doublet-TRD" ~ "Doublet-TRD",
                                  !is.na(selected_contig_TRB) &
                                    locus_TRB == "TRB" & 
                                    real_prod_TRB ~ "TRB",
                                  !is.na(selected_contig_TRD) &
                                    locus_TRD == "TRD" &
                                    real_prod_TRD ~ "TRD",
                                  .default = "No"),
          Light_Chain = case_when(selected_contig_TRA == "Doublet-TRA" ~ "Doublet-TRA",
                                  selected_contig_TRG == "Doublet-TRG" ~ "Doublet-TRG",
                                  !is.na(selected_contig_TRA) &
                                    locus_TRA == "TRA" & 
                                    real_prod_TRA ~ "TRA",
                                  !is.na(selected_contig_TRG) &
                                    locus_TRG == "TRG" &
                                    real_prod_TRG ~ "TRG",
                                  .default = "No"),
          Paired = case_when(Single_HC != "No" & Light_Chain != "No" ~ 1,
                              .default = 0)
    ) %>%
    relocate(Paired, Single_HC, Light_Chain, any_of(c("const_TRB", "const_TRA", "const_TRD", "const_TRG")), .after = cell_id) %>%
    relocate(any_of(c("v_mutation_pct_TRB", "v_mutation_pct_TRA", "v_mutation_pct_TRD", "v_mutation_pct_TRG",
            "cdr3_aa_length_TRB", "cdr3_aa_length_TRA", "cdr3_aa_length_TRD", "cdr3_aa_length_TRG")), .after = c_call_TRA) %>%
    # If this information is of interest, comment line
    select(-any_of(c("selected_contig_TRB", "selected_contig_TRA", "selected_contig_TRD", "selected_contig_TRG",
                     "locus_TRB", "locus_TRA", "locus_TRD", "locus_TRG",
                     "sequence_id_TRB", "sequence_id_TRA", "sequence_id_TRD", "sequence_id_TRG",
                     "real_prod_TRB", "real_prod_TRA", "real_prod_TRD", "real_prod_TRG")))
  write.table(anno.pair.wide, paste0("Processing_QC/", project_name,"_sadie_wide.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  paired <- anno.pair.wide$cell_id[anno.pair.wide$Paired > 0]

  # Filter unpaired (non-productive, liable and artifacts) from airr for changeo-10x.sh (without AssignGenes.py)(modified 2023.12.08)
  anno.pair %>%
    filter(cell_id %in% paired, real_prod) %>%
    select(-selected_contig, -real_prod) %>%
    write.table(paste0("Processing_QC/", project_name,"_ig.sadie_paired_airr.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

}
