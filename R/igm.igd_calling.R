#' IgM+IgD+ Calling
#'
#' Assign constant region as IgM.IgD to cells that have the same V,J and CDR3 (aa) instead of possibly calling them as Doublets in the pairing analysis.
#' @param singlets_prod_m Data frame in airr format.
#' @keywords pairing
#' @export
#'

igm.igd_calling <- function(singlets_prod_m, project_name){
  
  print("Contig selection with IgM+IgD+ option")

      # For IgM+IgD+ cells
  # Select cells with IgD contigs
  cells.dpos <- singlets_prod_m %>%
    group_by(cell_id) %>%
    filter(const == "IGHD")
  cells.dpos <- unique(cells.dpos$cell_id)
  # Select cells with IgM from IgD cells
  cells.dpos <- singlets_prod_m %>%
    group_by(cell_id) %>%
    filter(cell_id %in% cells.dpos, const == "IGHM")
  cells.dpos <- unique(cells.dpos$cell_id)

  # Check if v, j and cdr3 aa are the same
  cells.dpos <- singlets_prod_m %>%
    filter(cell_id %in% cells.dpos, const %in% c("IGHM", "IGHD")) %>%
    group_by(cell_id) %>%
    select(cell_id, v_call, j_call, cdr3_aa, const, c_call, sequence_id) %>%
    mutate(v_gene = str_split_i(v_call, "\\*", 1),
           same.v = ifelse(n_distinct(v_gene) == 1, TRUE, FALSE),
           same.j = ifelse(n_distinct(j_call) == 1, TRUE, FALSE)) %>%
    filter(same.v, same.j) %>%
    select(cell_id, cdr3_aa, const, sequence_id)

    # If zero cells are double positive, then return original input
    if(nrow(cells.dpos) > 0){
      # To wide
      cells.dpos <- cells.dpos %>%
        pivot_wider(id_cols = cell_id,
                    names_from = const,
                    names_glue = "{.value}_{const}",
                    values_from = c(cdr3_aa, sequence_id)) %>%
        mutate(same.cdr = ifelse(stringdist::stringdist(cdr3_aa_IGHM, cdr3_aa_IGHD, method = "lv") <= 3, TRUE, FALSE))

    # Save histogram to check if cutoff = 3 still valid
    pdf(paste0("Processing_QC/", project_name, "_HCDR3_Lv_dist_IgMIgD.pdf"))
    hist(stringdist::stringdist(cells.dpos$cdr3_aa_IGHM, cells.dpos$cdr3_aa_IGHD, method = "lv"), breaks = 20,
        main = project_name, xlab = "Levenshtein distance between HCDR3")
    dev.off()

    # Filter same cdr3 aa sequence
    cells.dpos <- cells.dpos %>%
      filter(same.cdr)

    # It's the same to select IgM or IgD since both are the same sequence, only C changes (we care more about M than D?)
    sequence_id.dpos <- c(cells.dpos$sequence_id_IGHM, cells.dpos$sequence_id_IGHD)
    selected.dpos <- rep(cells.dpos$sequence_id_IGHM, each = 2)
    cell_id.dpos <- cells.dpos$cell_id

    # Return to main df
    singlets_prod_m[singlets_prod_m$sequence_id %in% sequence_id.dpos, "const"] <- "IgM.IgD"
    singlets_prod_m[singlets_prod_m$sequence_id %in% sequence_id.dpos, "selected_contig"] <- selected.dpos
  }

  print(paste("Number of IgM+IgD+:", nrow(cells.dpos)))

  singlets_prod_m

}