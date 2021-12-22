#!/usr/bin/env Rscript
library(devtools)
library(phylostratr)
library(reshape2)
library(taxizedb)
library(dplyr)
library(readr)
library(magrittr)
library(ggtree)
library(knitr)
args <- commandArgs()
dir <- args[6]
setwd(dir)
protein_file <- args[7]
weights=uniprot_weight_by_ref()
focal_taxid <- '3702'
strata <-
  uniprot_strata(focal_taxid, from=2) %>%
  strata_apply(f=diverse_subtree, n=5, weights=weights) %>%
  use_recommended_prokaryotes %>%
  add_taxa(c('4932', '9606')) %>%
  uniprot_fill_strata
strata@data$faa[['3702']] <- protein_file