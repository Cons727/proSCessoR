#' Select Annotations from SADIE
#'
#' Select the annotations from SADIE runs and mark as liable the contigs that do not have complete sequences.
#' @param sadie_dir Path to directory of IgSADIE results. Default: Processing_QC/project.name_IgBLAST_SADIE_Runs
#' @param project_name Name of the project.
#' @keywords pairing
#' @export
#'


library(data.table)

IgSADIE <- function(sadie_dir = NULL, project_name) {
  
  if(is.null(sadie_dir)){
    sadie_dir <- paste0("Processing_QC/", project_name, "_IgBLAST_SADIE_Runs")
  }
  
  # 1. Get file names
  # From SADIE V ranges from -1 to -3, J ranges from -1 to -2, D is fixed at -1.
  # https://github.com/jwillis0720/sadie/blob/967cd1ddfbd6dfc5c67220839c52067211f596ae/src/sadie/airr/airr.py#L555
  # Liable from SADIE
  # https://github.com/jwillis0720/sadie/blob/967cd1ddfbd6dfc5c67220839c52067211f596ae/src/sadie/airr/airrtable/airrtable.py#L650C43-L650C43
  dt_grid <- data.table::data.table(
    penalty_rank      = 1:7,
    v_penalty         = c(-1, -1, -2, -2, -3, -3, -1),
    d_penalty         = c(-1, -1, -1, -1, -1, -1, -4),
    j_penalty         = c(-1, -2, -1, -2, -1, -2, -3),
    allow_vdj_overlap = c( F,  F,  F,  F,  F,  F,  T)
  )
  
  dt_grid[, penalties := sprintf("V%s_D%s_J%s", -v_penalty, -d_penalty, -j_penalty)]
  dt_grid[allow_vdj_overlap == TRUE, penalties := paste0(penalties, "_overlap")]
  dt_grid[, file := file.path(sadie_dir, paste0("SADIE_", penalties, ".tsv"))]
  
  # 2. Read all TSVs into a single list
  file_list <- lapply(dt_grid$file, data.table::fread)
  names(file_list) <- dt_grid$penalty_rank
  
  # 3. Combine into one data.table
  dt_all <- data.table::rbindlist(file_list, idcol = "penalty_rank")
  dt_all[, penalty_rank := as.integer(penalty_rank)]
  
  # Join the penalty metadata back to the sequences
  dt_all <- dt_grid[dt_all, on = "penalty_rank"]
  
  n_unique_seqs <- data.table::uniqueN(dt_all$sequence_id)
  cat("Annotating", n_unique_seqs, "sequences across", nrow(dt_grid), "penalty configurations...\n")
  
  # 4. Vectorized Liability Check
  cols <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")
  m <- dt_all[, lapply(.SD, function(x) is.na(x) | x == ""), .SDcols = cols]
  
  num_missing <- rowSums(m)
  num_cols <- length(cols)
  
  dt_all[, liable := data.table::fcase(
    num_missing == num_cols, TRUE, # all regions (FWR1 - FWR4) are missing
    num_missing == (num_cols - 1), TRUE, # only one region is found, 6 are missing
    num_missing == 0, FALSE, # the sequence has all regions
    !m$cdr2 & !m$fwr3 & m$cdr3, TRUE, # CDR2 present, FWR3 present, missing CDR3
    (m$fwr1 | m$cdr1 | m$fwr2 | m$cdr2) & !m$fwr3 & !m$cdr3, FALSE,
    # at least one V region is missing (FWR1 - CDR2), but FWR3 - CDR3 is present
    !m$fwr3 & !m$cdr3 & m$fwr4, TRUE, # FWR3 and CDR3 are present, but FWR4 is missing
    default = TRUE
  )]
  
  # 5. If liable == FALSE, we prioritize the earliest penalty (e.g., rank 1).
  # If liable == TRUE (across all configs), the loop would retain the *last* penalty (e.g., rank 7).
  # fifelse creates a temporary custom ranking to force this behavior during the sort.
  max_rank <- nrow(dt_grid)
  dt_all[, sort_rank := data.table::fifelse(
    liable == FALSE, 
    penalty_rank, 
    (max_rank + 1L) - penalty_rank
  )]
  
  # Sort and extract the top "winning" row for each sequence
  data.table::setorder(dt_all, sequence_id, liable, sort_rank)
  anno <- unique(dt_all, by = "sequence_id")
  anno[, sort_rank := NULL]
  
  # 6. Final Outputs
  cat(sprintf("%s sequences still liable\n", nrow(anno[liable == TRUE])))
  print(table(anno$penalties))
  
  stopifnot(nrow(anno) == n_unique_seqs)
  
  return(anno)
}