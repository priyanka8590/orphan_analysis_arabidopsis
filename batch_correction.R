#if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("sva")
#library("sva")
#library("ggplot2")
#library("gridExtra")
library("devtools")
devtools::install_github("zhangyuqing/sva-devel")
setwd("~/Box/analysis_orphan_gene_aim_3/")
count_matrix <- read.table("final_expression_table_with_high_depth_25.tsv", sep = "\t", header = TRUE)
count_matrix <- as.data.frame(count_matrix)
#trans_count <- t(count_matrix)
studies_with_one_run <- read.delim("runs_to_be_removed")
studies_with_one_run <- as.data.frame(studies_with_one_run)
new_variable <- as.vector(studies_with_one_run$runs)
df <- count_matrix[,!(names(count_matrix) %in% new_variable)]
#count_matrix <- count_matrix[!(count_matrix$rowname %in% studies_with_one_run$study),]
#uncorrected_data=read.table("QNnorm.tsv", header=TRUE, sep="\t", as.is = c(1))
#uncorrected_data_matrix=as.matrix(uncorrected_data)
#metadata_table=read.table("final_metadata_for_1301_runs", header = TRUE, sep = "\t", fill = TRUE)
read.delim("~/Box Sync/analysis_orphan_gene_aim_3/testing_part_two/final_metadata_for_1301_runs")
metadata_table <- read.delim("final_metadata_for_4741_runs_for_deseq2_original.txt")
metadata_table <- metadata_table[!(metadata_table$run_accession %in% new_variable),]
batches<-as.vector(metadata_table$study_accession)
#condition<-metadata_table$condition
#group<-sapply(as.character(condition), switch, "aerial" = 1, "non-aerial" = 2, USE.NAMES = F)
df_mat <- data.matrix(df)
row.names(df_mat) <- df$rowname
df_mat <- df_mat[,c(-1)]
count_matrix_mat <- data.matrix(count_matrix)
final_matrix <- cbind(corrected_data, count_matrix[,names(count_matrix) %in% new_variable])
corrected_data<-ComBat_seq(counts=df_mat, batch=batches)

#DESeq2 normalization
library(DESeq2)

