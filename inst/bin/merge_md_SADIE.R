#!/usr/bin/env Rscript

library(tidyverse)

project_name <- commandArgs(trailingOnly = TRUE)[1]

# Merge ChangeO with IgSADIE
# From merge_ChangeO()
anno.pair.wide <- read.table(paste0("Processing_QC/", project_name,"_sadie_wide.tsv"),
                             header = TRUE, sep = "\t", fill = TRUE, check.names = FALSE)

md <- read.csv(paste0("Processing_QC/", project_name, "_metadata_with_bc.csv"))
# Use lubridate::today() to add date
out.name <- paste0(lubridate::today(),"_", project_name, "_IgSADIE.tsv")

left_join(anno.pair.wide, md, by = join_by(cell_id == barcode)) %>%
  relocate(Batch_Pool, Sort_Population, PTID, Visit, Status, HTO_Assignment, .after = cell_id) %>%
  write.table(out.name, sep = "\t", row.names = F, quote = F)