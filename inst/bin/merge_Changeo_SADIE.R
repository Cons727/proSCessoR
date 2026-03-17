#!/usr/bin/env Rscript

library(tidyverse)

project_name <- commandArgs(trailingOnly = TRUE)[1]

# Merge ChangeO with IgSADIE
# From merge_ChangeO()
anno.pair.wide <- read.table(paste0("Processing_QC/", project_name,"_sadie_wide.tsv"),
                             header = TRUE, sep = "\t", fill = TRUE, check.names = FALSE)
changeo <- read.table(paste0("Processing_QC/", project_name, "_changeo/", project_name, "_heavy_clone-pass.tsv"),
                      header = TRUE, sep = "\t", fill = TRUE, check.names = FALSE)

ig.changeo <- anno.pair.wide %>%
  left_join(select(.data = changeo, cell_id, clone_id), by = "cell_id") %>%
  group_by(clone_id) %>%
  mutate(clone_size = n()) %>%
  ungroup() %>%
  mutate(clone_size = ifelse(is.na(clone_id), NA, clone_size)) %>%
  relocate(clone_id, clone_size, .after = Light_Chain)

md <- read.csv(paste0("Processing_QC/", project_name, "_metadata_with_bc.csv"))
# Use lubridate::today() to add date
out.name <- paste0(lubridate::today(),"_", project_name, "_ChangeO-IgSADIE.tsv")

left_join(ig.changeo, md, by = join_by(cell_id == barcode)) %>%
  relocate(Batch_Pool, Sort_Population, PTID, Visit, Status, HTO_Assignment, .after = cell_id) %>%
  write.table(out.name, sep = "\t", row.names = F, quote = F)