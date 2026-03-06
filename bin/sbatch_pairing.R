#!/usr/bin/env Rscript

devtools::load_all("~/proSCessoR")
library(tidyverse)

project.name <- commandArgs(trailingOnly = TRUE)[1]
vdj.type <- commandArgs(trailingOnly = TRUE)[2]

batches <- grep("^P[[:digit:]]", list.dirs(path = ".", full.names = FALSE, recursive = FALSE), value = TRUE)

anno <- IgSADIE(project_name = project.name)

if(vdj.type == "vdj_b"){
    pairing(anno = anno, project_name = project.name, kl.light = FALSE)
}
if(vdj.type == "vdj_t"){
    tcr.pairing(anno = anno, project_name = project.name)
}