file_to_be_normalized <- read.table(file = "deseq2norm_25.tsv", sep = "\t", header = TRUE)
file_to_be_normalized_df <- as.data.frame(file_to_be_normalized)
#NumReads_matSCN <- round(file_to_be_normalized_df[,-1])
NumReads_matSCN <- file_to_be_normalized_df[,-1]
NumReads_matSCN <- round(file_to_be_normalized_df[,-1])
rownames(NumReads_matSCN) <- file_to_be_normalized_df[,1]
NumReads_matSCN <- data.matrix(NumReads_matSCN)
metadata_table <- read.delim("~/Box/analysis_orphan_gene_aim_3/final_metadata_for_4741_runs_for_deseq2_original.txt")
NumReads_matSCN <- NumReads_matSCN + 1
condition <- metadata_table$condition
coldata<-data.frame(row.names=colnames(NumReads_matSCN),condition)
dds <- DESeqDataSetFromMatrix(countData = round(NumReads_matSCN), colData = coldata,design=~condition)
dds<-estimateSizeFactors(dds) 
norm_deseq2 <- counts(dds,normalized=TRUE)# Normalized counts
norm_deseq2 <- as.data.frame(norm_deseq2)
rownames(norm_deseq2) <- rownames(NumReads_matSCN)
norm_deseq2 <- tibble::rownames_to_column(norm_deseq2)
colnames(norm_deseq2) <- colnames(file_to_be_normalized)
write.table(norm_deseq2, file = "deseq2norm_25.tsv", sep = "\t", quote = F, row.names = F)