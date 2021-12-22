setwd("/work/LAS/mash-lab/bhandary/analysis_orphan_genes_in_arabidopsis")
#setwd("/Users/bhandary/Box\ Sync/analysis_orphan_genes_in_arabidopsis")
#library("edgeR")
#library("dplyr")
#library("tibble")
#library("ggplot2")
#devtools::install_github("rhondabacher/SCnorm")
library("SCnorm")
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#library('DESeq2')

############################################################################################
# Setup variables
############################################################################################

args <- commandArgs()
input_filename<-args[6]
output_filename<-args[7]
file_to_be_normalized <- read.table(input_filename, sep = "\t", header = TRUE)
file_to_be_normalized_df <- as.data.frame(file_to_be_normalized)
#file_to_be_normalized <- read.table(file = "final_counts_to_check_best_norm", sep = "\t", header = TRUE)
#SC normalization
NumReads_matSCN <- round(file_to_be_normalized_df[,-1])
rownames(NumReads_matSCN) <- file_to_be_normalized_df[,1]
NumReads_matSCN <- data.matrix(NumReads_matSCN)
NumReads_matSCN[NumReads_matSCN > 0] <- floor(NumReads_matSCN[NumReads_matSCN > 0]) #first remove genes that have no counts in any samples
NumReads_matSCN_final <- NumReads_matSCN[rowSums(NumReads_matSCN) > 10,] #Only retain rows that have a positive total expression per gene
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