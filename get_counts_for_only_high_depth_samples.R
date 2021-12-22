setwd("~/Box Sync/analysis_orphan_gene_aim_3")
exp_data <- read.table(file="salmon_consolidated_all_5210.tsv", sep = "\t", header=TRUE)
transposed_exp_data <- t(exp_data)
high_depth_samples <- read.table(file="low_dep")