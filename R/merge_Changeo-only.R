#' Merge Change-O Output with enclone-VQuest
#'
#' Merge heavy_clone-pass.tsv with VQuest-enclone data frame. File should not contain barcode overlaps. Returns a list of data frames and saves output as tsv and RData.
#' @param clone_path Change-O output path.
#' @param v_enclone_path Path to v-enclone.tsv
#' @param md_path Path to project.name_metadata_with_bc.csv; must include barcode. Default to proSCessoR's output.
#' @param project_name Project name.
#' @keywords changeo
#' @export
#'

merge_Changeo <- function(clone_path, v_enclone_path, project_name, md_path = NULL){

  changeo <- read.table(clone_path, header = TRUE, sep = "\t", fill = TRUE, check.names = FALSE)
  v_enclone <- read.table(v_enclone_path, header = TRUE, sep = "\t", fill = TRUE, check.names = FALSE)

  c_v_e <- v_enclone %>%
    left_join(select(changeo, cell_id, clone_id), by = c("barcode" = "cell_id"))
  c_v_e <- c_v_e %>%
    group_by(clone_id) %>%
    mutate(changeo_clone_ncells = n()) %>%
    ungroup() %>%
    mutate(changeo_clone_ncells = ifelse(is.na(clone_id), NA, changeo_clone_ncells)) %>%
    relocate(clone_id, .after = Light_Chain) %>%
    relocate(changeo_clone_ncells, .after = clone_id)

  if(is.null(md_path)){
    md <- read.csv(paste0("Processing_QC/", project_name, "_metadata_with_bc.csv"))
  } else{
    md <- read.csv(md_path)
  }
  
  c_v_e <- md %>%
    left_join(c_v_e, ., by = c(barcode = "barcode")) %>%
    relocate(c(Batch_Pool, Sort_Population, PTID, Visit, Status, HTO_Assignment), .after = barcode) %>%
    relocate(Paired, .before = Unpaired) %>%
    mutate(Top_contig_HC = ifelse(Top_contig_HC == 0, NA, Top_contig_HC)) %>%
    mutate(Top_contig_KC = ifelse(Top_contig_KC == 0, NA, Top_contig_KC)) %>%
    mutate(Top_contig_LC = ifelse(Top_contig_LC == 0, NA, Top_contig_LC))
  # Use lubridate::today() to add date
  out.name <- paste0(lubridate::today(),"_", project_name, "_ChangeO-Vquest.tsv")
  write.table(c_v_e, out.name, sep = "\t", row.names = F, quote = F)

  c_v_e
}
