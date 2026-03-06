#' Contig Selection
#'
#' Determines best contigs for each cell using method. Saves fasta and annotation files for productive singlets.
#' @param batches Vector of batch or pool names.
#' @param project_name Name of the project.
#' @param method Heuristic to select best contigs. Can be "top3" for first 3 highest contigs based on number of reads, or "top2" for first 2 highest contigs based on reads for each locus. Default: "top3".
#' @param file IgBLAST output file path in airr format. Default: Processing_QC/singlets_contigs/project.name_singlets_contigs_igblast.tsv
#' @param igm.igd Logical. Call IgM+IgD+ instead of Doublet. Default FALSE.
#' @param vdj_type String. VDJ type, either "vdj_b" for BCR or "vdj_t" for TCR. Default "vdj_b".
#' @keywords pairing
#' @export
#'

contig_selection <- function(batches, project_name, igm.igd = FALSE, vdj_type = "vdj_b", method = "top3", file = NULL){

  if(is.null(file)){
    file <- paste0("Processing_QC/singlets_contigs/", project_name, "_singlets_contigs_igblast.tsv")
  }

  singlets_prod <- read.delim(file, sep = "\t") %>%
    filter(productive)
  
  singlets_prod$const <- stringr::str_split_i(singlets_prod$c_call, "\\*", 1)
  singlets_prod$cell_id <- stringr::str_split_i(singlets_prod$sequence_id, "_contig_", 1)

  # Add read and umi counts info from all_contig_annotation.csv
  cts <- paste0(batches, "/multi/", vdj_type)
  cts <- sapply(cts, get_CR_file, file = "all_contig_annotations.csv")
  cts <- lapply(cts, read.csv)
  names(cts) <- batches
  cts <- lapply(seq_along(cts),
                function(i) cts[[i]] %>%
                  select(contig_id, reads, umis) %>%
                  mutate(contig_id = paste0(names(cts[i]), "_", contig_id))
  )
  cts <- do.call(rbind, cts)

  singlets_prod <- left_join(singlets_prod, cts, by = join_by(sequence_id == contig_id)) %>%
    relocate(cell_id, .before = sequence_id)
  rm(cts)

  # Top 3 contigs per cell
  if(method == "top3"){
    singlets_prod_m <- singlets_prod %>%
      group_by(cell_id) %>%
      arrange(desc(reads), .by_group = TRUE) %>%
      slice_max(n = 3, order_by = reads, with_ties = FALSE) %>%
      group_by(cell_id, locus) %>%
      mutate(ratio = min(reads) / max(reads),
            selected_contig = ifelse(min(reads) / max(reads) < 0.5 | n() == 1, sequence_id[which.max(reads)], paste0("Doublet-", locus))) %>%
      ungroup()
  }

  # Top 2 contigs for each locus per cell
  if(method == "top2"){
    singlets_prod_m <- singlets_prod %>%
      group_by(cell_id, locus) %>%
      arrange(desc(reads), .by_group = TRUE) %>%
      slice_max(order_by = reads, n = 2) %>%
      mutate(selected_contig = ifelse(min(reads) / max(reads) < 0.5 | n() == 1, sequence_id[which.max(reads)], paste0("Doublet-", locus))) %>%
      ungroup()
  }

  if(igm.igd){
    singlets_prod_m <- igm.igd_calling(singlets_prod_m = singlets_prod_m, project_name = project_name)
  }

  # It doesn't matter where the doublet is coming from, it will be called as Doublet
  doublets.cells <- unique(singlets_prod_m$cell_id[grepl("Doublet", singlets_prod_m$selected_contig)])
  doublets.df <- singlets_prod_m %>%
    filter(cell_id %in% doublets.cells, # keep doublets cells
          grepl("Doublet", selected_contig)) %>% # keep doublet rows
    distinct(cell_id, .keep_all = TRUE) # keep highest doublet row

  # Since dataframe is sorted by highest reads, first unique values will be highest
  singlets_prod_f <- singlets_prod_m %>%
    filter(!cell_id %in% doublets.cells) %>%
    distinct(cell_id, selected_contig, .keep_all = TRUE) %>%
    rbind(., doublets.df) %>%
    mutate(sequence_id = paste0(cell_id, "_", locus))
  write.table(singlets_prod_f, paste0("Processing_QC/", project_name, "_singlets_prod_contig.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  # Create fasta of singlets
  fasta.string <- singlets_prod_f$sequence[!singlets_prod_f$cell_id %in% doublets.cells]
  names(fasta.string) <- singlets_prod_f$sequence_id[!singlets_prod_f$cell_id %in% doublets.cells]

  fasta.string <- Biostrings::DNAStringSet(fasta.string)
  Biostrings::writeXStringSet(fasta.string, paste0("Processing_QC/", project_name, "_singlets_prod_contig.fasta"))

}