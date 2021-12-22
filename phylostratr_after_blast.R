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
num_threads <- args[8]
weights=uniprot_weight_by_ref()
focal_taxid <- '3702'
strata <-
  uniprot_strata(focal_taxid, from=2) %>%
  strata_apply(f=diverse_subtree, n=5, weights=weights) %>%
  use_recommended_prokaryotes %>%
  add_taxa(c('4932', '9606')) %>%
  uniprot_fill_strata
strata@data$faa[['3702']] <- protein_file
strata %>% strata_convert(target='all', to='name') %>% sort_strata %>% plot
strata <- strata_blast(strata, blast_args=list(nthreads=num_threads)) %>% strata_besthits
results <- merge_besthits(strata)
phylostrata <- stratify(results, classify_by_adjusted_pvalue(0.001))
write.csv(phylostrata, file=paste0(protein_file,"_phylostrata_table.csv"))