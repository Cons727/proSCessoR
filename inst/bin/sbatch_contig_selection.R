#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(proSCessoR))
suppressPackageStartupMessages(library(tidyverse))

project.name <- commandArgs(trailingOnly = TRUE)[1]
vdj.type <- commandArgs(trailingOnly = TRUE)[2]
igm.igd <- commandArgs(trailingOnly = TRUE)[3]

batches <- grep("^P[[:digit:]]", list.dirs(path = ".", full.names = FALSE, recursive = FALSE), value = TRUE)

contig_selection(batches = batches, project_name = project.name, vdj_type = vdj.type, igm.igd = igm.igd, method = "top3")
