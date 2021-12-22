#setwd("/work/LAS/mash-lab/bhandary/analysis_orphan_genes_in_arabidopsis/testing_part_two")
setwd("~/Box/analysis_orphan_gene_aim_3/testing_part_two")
setwd("~/Box/analysis_orphan_gene_aim_3/")
library("edgeR")
#library("dplyr")
#library("tibble")
#library("ggplot2")
#devtools::install_github("rhondabacher/SCnorm")
library("SCnorm")
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library('DESeq2')


library('preprocessCore')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("preprocessCore")

############################################################################################
# Setup variables
############################################################################################

args <- commandArgs()
input_filename<-args[6]
output_filename<-args[7]
file_to_be_normalized <- read.table(input_filename, sep = "\t", header = TRUE)
file_to_be_normalized_df <- as.data.frame(file_to_be_normalized)
file_to_be_normalized <- read.table(file = "final_counts_for_testing_part_two_tab.txt", sep = "\t", header = TRUE)
file_to_be_normalized <- read.table(file = "final_expression_table_with_high_depth.tsv", sep = "\t", header = TRUE)
file_to_be_normalized <- read.table(file = "expression_table_with_fap_samples.tsv", sep = "\t", header = TRUE)
file_to_be_normalized <- read.table(file = "final_expression_table_with_high_depth_25.tsv", sep = "\t", header = TRUE)
file_to_be_normalized <- read.table(file = "deseq2norm_25.tsv", sep = "\t", header = TRUE)
#SC normalization
#final_expression_table_with_high_depth.tsv
file_to_be_normalized_df <- as.data.frame(file_to_be_normalized)
#NumReads_matSCN <- round(file_to_be_normalized_df[,-1])
NumReads_matSCN <- file_to_be_normalized_df[,-1]
rownames(NumReads_matSCN) <- file_to_be_normalized_df[,1]
NumReads_matSCN <- data.matrix(NumReads_matSCN)
NumReads_matSCN[NumReads_matSCN > 0] <- floor(NumReads_matSCN[NumReads_matSCN > 0]) #first remove genes that have no counts in any samples
NumReads_matSCN <- NumReads_matSCN[rowSums(NumReads_matSCN) > 1,] #Only retain rows that have a positive total expression per gene
#Conditions <- factor(c(rep("S",ncol(file_to_be_normalized_df)/2),rep("N",ncol(file_to_be_normalized_df)/2)))
#print(Conditions)
#Conditions <- factor(c(rep("S",80),rep("N",80)))
#Conditions <- c(rep(1, 35), rep(2, 6), rep(3, 12), rep(4,1), rep(5, 8), rep(6, 3), rep(7, 5), rep(8, 38), rep(9, 3), rep(10, 3), rep(11, 45))
#Conditions <- c(rep(1, 35), rep(2, 6), rep(3, 12), rep(4, 8), rep(5, 3), rep(6, 5), rep(7, 38), rep(8, 3), rep(9, 3), rep(10, 45))
Conditions <- c(rep(1, 158))
#NumReads_SCNorm <- SCnorm(Data = NumReads_matSCN, Conditions = Conditions, PrintProgressPlots = TRUE, NCores = 3, useZerosToScale = TRUE, FilterCellNum = 15, FilterExpression = 0)
countDeptEst <- plotCountDepth(Data = NumReads_matSCN_final, Conditions = Conditions, FilterCellProportion = 0.1, NCores = 3, FilterExpression=3)
NumReads_SCNorm <- SCnorm(Data = NumReads_matSCN_final, Conditions = Conditions, PrintProgressPlots = TRUE, NCores = 3, PropToUse = 0.1, FilterCellNum = 10, ditherCounts = TRUE, Thresh = 0.1)
write.table(NumReads_SCNorm, file = "SCnorm.tsv", sep = "\t", quote = F)

#Quantile normalization
NumReads_matSCN <- round(file_to_be_normalized_df[,-1])
rownames(NumReads_matSCN) <- file_to_be_normalized_df[,1]
NumReads_matSCN <- data.matrix(NumReads_matSCN)
NumReads_mat_qn_norm <- normalize.quantiles(NumReads_matSCN, copy = TRUE)
NumReads_mat_qn_norm <- as.data.frame(NumReads_mat_qn_norm)
rownames(NumReads_mat_qn_norm) <- rownames(NumReads_matSCN)
NumReads_mat_qn_norm <- tibble::rownames_to_column(NumReads_mat_qn_norm)
colnames(NumReads_mat_qn_norm) <- colnames(file_to_be_normalized)
write.table(NumReads_mat_qn_norm, file = "QNnorm.tsv", sep = "\t", quote = F, row.names = F)
                                   
#DESeq2
NumReads_matSCN <- round(file_to_be_normalized_df[,-1])
rownames(NumReads_matSCN) <- file_to_be_normalized_df[,1]
NumReads_matSCN <- data.matrix(NumReads_matSCN)
#metadata_table <- read.delim("~/Box/analysis_orphan_gene_aim_3/testing_part_two/final_metadata_for_1301_runs")
metadata_table <- read.delim("~/Box/analysis_orphan_gene_aim_3/final_metadata_for_4741_runs_for_deseq2.txt")
metadata_table <- read.delim("~/Box/analysis_orphan_gene_aim_3/complete_run_metadata_athaliana_pb_4786runs_08242021.txt")
metadata_table <- read.delim("~/Box/analysis_orphan_gene_aim_3/final_metadata_for_4741_runs_for_deseq2_original.txt")
#coldata <- colnames(NumReads[2:2955]) 
#condition<-factor(c(rep("S",99)))
#condition<-factor(c(rep("S",1477),rep("N",1477)))
##condition <- c(rep(1, 35), rep(2, 6), rep(3, 12), rep(4, 8), rep(5, 3), rep(6, 5), rep(7, 38), rep(8, 3), rep(9, 3), rep(10, 45)) #for first testing
condition <- as.factor(c(rep(1,8),rep(2,6),rep(3,36),rep(4,17),rep(5,1),rep(6,234),rep(7,18),rep(8,82),rep(9,58),rep(10,18),rep(11,636),rep(12,20),rep(13,6),rep(14,96),rep(15,2),rep(16,9),rep(17,4),rep(18,9),rep(19,31),rep(20,6),rep(21,6)))
#condition <- c(rep(1, 158))
#condition <- factor(c(rep("S",79),rep("N",79)))
NumReads_matSCN <- NumReads_matSCN + 1
condition <- metadata_table$condition
batches=as.vector(metadata_table$study_accession)
coldata<-data.frame(row.names=colnames(NumReads_matSCN),condition,batches) 
coldata<-data.frame(row.names=colnames(NumReads_matSCN),condition)
dds <- DESeqDataSetFromMatrix(countData = round(NumReads_matSCN), colData = coldata,design=~condition)
dds<-estimateSizeFactors(dds) 
norm_deseq2 <- counts(dds,normalized=TRUE)# Normalized counts
norm_deseq2 <- as.data.frame(norm_deseq2)
rownames(norm_deseq2) <- rownames(NumReads_matSCN)
norm_deseq2 <- tibble::rownames_to_column(norm_deseq2)
colnames(norm_deseq2) <- colnames(file_to_be_normalized)
write.table(norm_deseq2, file = "deseq2norm_25.tsv", sep = "\t", quote = F, row.names = F)

#TMM
NumReads_matSCN <- round(file_to_be_normalized_df[,-1])
rownames(NumReads_matSCN) <- file_to_be_normalized_df[,1]
NumReads_matSCN <- data.matrix(NumReads_matSCN)
#count_data$Gene_ID<-NULL
#dim(count_data)
#group<-rep("dummy",dim(count_data)[2])
##group <- c(rep(1, 35), rep(2, 6), rep(3, 12), rep(4, 8), rep(5, 3), rep(6, 5), rep(7, 38), rep(8, 3), rep(9, 3), rep(10, 45))  #for first testing
group <- as.factor(c(rep(1,8),rep(2,6),rep(3,36),rep(4,17),rep(5,1),rep(6,234),rep(7,18),rep(8,82),rep(9,58),rep(10,18),rep(11,636),rep(12,20),rep(13,6),rep(14,96),rep(15,2),rep(16,9),rep(17,4),rep(18,9),rep(19,31),rep(20,6),rep(21,6)))
#import raw counts into EdgeR
d <- DGEList(counts=NumReads_matSCN, group=group)
#may have runs without any expression, just discard them
keep <- d$samples$lib.size != 0
d <- d[ ,keep]
no <- d$samples$lib.size == 0
colnames(d[ ,no])
TMM <- calcNormFactors(d, method="TMM")
CPM<-cpm(TMM, normalized.lib.sizes = T)
logcpm <- log(CPM+1)
norm_counts.table.TMM <- t(t(d$counts)*(TMM$samples$norm.factors))
norm_counts.table.TMM <- as.data.frame(norm_counts.table.TMM)
rownames(norm_counts.table.TMM) <- rownames(NumReads_matSCN)
norm_counts.table.TMM <- tibble::rownames_to_column(norm_counts.table.TMM)
colnames(norm_counts.table.TMM) <- colnames(file_to_be_normalized)
write.table(norm_counts.table.TMM, file="./normalizedCountsTMM.txt", sep="\t", quote=F, row.names = F)
write.table(logcpm, file = "logcpm.tsv", sep = "\t", quote = F)

#Smooth Quantile
#install.packages("qsmooth")
#library(quantro)
library("qsmooth")
NumReads_matSCN <- round(file_to_be_normalized_df[,-1])
rownames(NumReads_matSCN) <- file_to_be_normalized_df[,1]
NumReads_matSCN <- data.matrix(NumReads_matSCN)
#group <- c(rep(1, 35), rep(2, 6), rep(3, 12), rep(4, 8), rep(5, 3), rep(6, 5), rep(7, 38), rep(8, 3), rep(9, 3), rep(10, 45))  #for first testing
group <- as.factor(c(rep(1,8),rep(2,6),rep(3,36),rep(4,17),rep(5,1),rep(6,234),rep(7,18),rep(8,82),rep(9,58),rep(10,18),rep(11,636),rep(12,20),rep(13,6),rep(14,96),rep(15,2),rep(16,9),rep(17,4),rep(18,9),rep(19,31),rep(20,6),rep(21,6)))
qsmooth_object <- qsmooth(object = NumReads_matSCN, group_factor = group)
matdensity(log2(qsmoothData(qsmooth_object)+1), groupFactor = group,
           main = "qsmooth normalized data",
           xlab = "Expression (log2 scale)", ylab = "density")
legend('topright', levels(factor(group)), col = 1:2, lty = 1)

norm_qsmooth <- as.data.frame(qsmoothData(qsmooth_object))
rownames(norm_qsmooth) <- rownames(NumReads_matSCN)
norm_qsmooth <- tibble::rownames_to_column(norm_qsmooth)
colnames(norm_qsmooth) <- colnames(file_to_be_normalized)
write.table(norm_qsmooth, file = "qsmooth_norm.tsv", sep = "\t", quote = F, row.names = F)
