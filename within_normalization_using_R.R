#setwd("/Users/bhandary/Box Sync/analysis_orphan_gene_aim_3")
setwd("/work/LAS/mash-lab/bhandary/batch_normalization/pbs_scripts")
library("edgeR")
library("dplyr")
library("stringr")

############################################################################################
# Setup variables
############################################################################################

args <- commandArgs()
input_filename<-args[6]
modified_input_filename = str_remove(input_filename, "_for_norm_within_studies.tsv")
#output_filename<-args[7]
#testing <- read.table("arabidopsis_gene_attributes.txt", sep = "\t", header = TRUE)
#testing_df <- as.data.frame(testing)
#testing_df$length <- abs(testing_df$Gene.start..bp. - testing_df$Gene.end..bp.)
#file_to_be_normalized <- read.table(file = "DRP005813_counts_for_norm_within_studies.tsv", sep = "\t", header = TRUE)
file_to_be_normalized <- read.table(input_filename, sep = "\t", header = TRUE)
file_to_be_normalized_df <- as.data.frame(file_to_be_normalized)
NumReads_matSCN <- round(file_to_be_normalized_df[,-1])
rownames(NumReads_matSCN) <- file_to_be_normalized_df[,1]
NumReads_matSCN <- data.matrix(NumReads_matSCN)

#Size normalization
tpm_norm <- as.data.frame(apply(NumReads_matSCN, 2, function(x) (x - min(x))/(max(x) - min(x))))
#tpm_norm$Gene <- file_to_be_normalized_df$Gene
#tpm_sizenorm <- tpm_norm %>% select(Gene, everything())
tpm_norm <- tibble::rownames_to_column(tpm_norm)
colnames(tpm_norm) <- colnames(file_to_be_normalized)
write.table(tpm_norm, file = paste0(modified_input_filename, "_sizenorm.tsv"), sep = "\t", quote = F, row.names = F)

#write.table(tpm_norm, file = "testing_sizenorm.tsv", sep = "\t", quote = F, row.names = F)

#CPM
group<-rep("dummy",dim(NumReads_matSCN)[2])
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
logcpm <- as.data.frame(logcpm)
logcpm <- tibble::rownames_to_column(logcpm)
colnames(logcpm) <- colnames(file_to_be_normalized)
write.table(logcpm, file=paste0(modified_input_filename, "_cpmnorm.tsv"), sep="\t", quote=F, row.names = F)
#write.table(logcpm, file = "testing_cpm.tsv", sep = "\t", quote = F, row.names = F)

#TPM
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

