library(DESeq2)
#install.packages("pasilla")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pasilla")
library("pasilla")


pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata",
                       "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE)


cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

setwd("~/Box/analysis_orphan_gene_aim_3/deseq_analysis")
setwd("~/Box/analysis_orphan_gene_aim_3")
salmon_counts <- read.table(file = "salmon_consolidated_all_5210.tsv", sep = "\t", header = TRUE)
count_matrix <- read.table("final_expression_table_with_high_depth_25.tsv", sep = "\t", header = TRUE)
studies_with_one_run <- read.delim("runs_to_be_removed")
studies_with_one_run <- as.data.frame(studies_with_one_run)
new_variable <- as.vector(studies_with_one_run$runs)
df <- count_matrix[,!(names(count_matrix) %in% new_variable)]
metadata_table <- read.delim("final_metadata_for_4741_runs_for_deseq2_original.txt")
metadata_table <- metadata_table[!(metadata_table$run_accession %in% new_variable),]

salmon_counts <- subset(salmon_counts, Gene %in% df$rowname) #to get filtered genes - 53052 genes , df is final_expression file dataframe
salmon_counts <- as.data.frame(salmon_counts)
salmon_counts_pi_def <- salmon_counts %>% select(Gene, SRR3087771, SRR3087772, SRR3087773, SRR3087774, SRR3087775, SRR3087776, SRR3087777, SRR3087778, SRR3087779, SRR3087780, SRR3087781, SRR3087782, SRR3087783, SRR3087784, SRR3087785, SRR3087786)
salmon_counts_pi_def_mat <- data.matrix(salmon_counts_pi_def)
rownames(salmon_counts_pi_def_mat) <- salmon_counts_pi_def[,1]
salmon_counts_pi_def_mat <- salmon_counts_pi_def_mat[,-1]
salmon_counts_pi_def_mat <- round(salmon_counts_pi_def_mat)

#To get metadata
library(dplyr)
#dea_anno <- filter(metadata_table, study_accession %in% 'SRP068145')
#dea_exp1_anno <- filter(metadata_table, study_accession %in% 'SRP082470')
dea_exp1_anno <- filter(metadata_table, study_accession %in% 'SRP149149')
drop <- c("treatment", "genotype")
dea_exp1_anno <- dea_exp1_anno[,!(names(dea_exp1_anno) %in% drop)]
#dea_exp1_anno$treatment <- c("Col_L_3", "Col_L_2", "Col_L_1", "Col_D_3", "Col_D_2", "Col_D_1", "rrc_L_3", "rrc_L_3", "rrc_L_2", "rrc_L_2", "rrc_L_1", "rrc_L_1", "rrc_D_3", "rrc_D_3", "rrc_D_2", "rrc_D_2", "rrc_D_1", "rrc_D_1")
dea_exp1_anno$treatment <- c("Light1", "Light2", "Light3", "Dark1", "Dark2", "Dark3", "Light1", "Light1", "Light2", "Light2", "Light3", "Light3", "Dark1", "Dark1", "Dark2", "Dark2", "Dark3", "Dark3")
dea_exp1_anno$genotype <- c("Col0", "Col0", "Col0", "Col0", "Col0", "Col0", "rrc1", "rrc1", "rrc1", "rrc1", "rrc1", "rrc1", "rrc1", "rrc1", "rrc1", "rrc1", "rrc1", "rrc1")
rownames(dea_exp1_anno) <- dea_exp1_anno[,1]
dea_exp1_anno <- dea_exp1_anno[,-1]
dea_exp1_anno <- dea_exp1_anno[order(row.names(dea_exp1_anno)),]
dea_exp1_anno <- dea_exp1_anno[,c("treatment", "genotype", "run_accession")]
#or alternatively,if metadata table doesn't have some runs
dea_exp1_anno <- read_excel("dea_exp1_anno.xlsx")
write.table(dea_exp1_anno, file = "metadata_for_dea_exp1", sep = "\t", quote = F, row.names = F)
dea_exp1_anno <- read.delim("metadata_for_dea_exp1", row.names = 1)
dea_exp1_anno <- dea_exp1_anno[order(row.names(dea_exp1_anno)),]
dea_exp1_anno <- dea_exp1_anno[,c("tissue", "treatment")]

dea_exp1_anno$genotype <- factor(dea_exp1_anno$genotype)
dea_exp1_anno$treatment <- factor(dea_exp1_anno$treatment)
dea_exp1_anno$run_accession <- factor(dea_exp1_anno$run_accession)
dds <- DESeqDataSetFromMatrix(countData = salmon_counts_pi_def_mat, colData = dea_exp1_anno, design = ~ treatment + genotype + treatment:genotype)
#keep <- rowSums(counts(dds)) >= 10
salmon_counts_pi_def_mat[salmon_counts_pi_def_mat > 0] <- floor(salmon_counts_pi_def_mat[salmon_counts_pi_def_mat > 0]) #first remove genes that have no counts in any samples
#salmon_counts_pi_def_mat <- salmon_counts_pi_def_mat[rowSums(salmon_counts_pi_def_mat) > 1,] #never use this
#dds$treatment <- factor(dds$treatment, levels = c("Phosphate deficiency", "Phosphate sufficiency"))
#dds$treatment <- factor(dds$treatment, levels = c("3hrs red light irradiation (7Ã«_molm-2s-1)", "none"))
#ddsColl <- collapseReplicates(dds, dds$treatment, dds$run_accession, renameCols = TRUE)
#dds$group <- factor(paste0(dds$genotype, dds$treatment))
#design(dds) <- ~ group
#dds$treatment <- factor(dds$treatment, levels = c("Col_D", "Col_L", "rrc_L", "rrc_D"))
dea_exp1_anno$ind.n <- factor(rep(rep(1:3, each = 3), 3))
dds <- DESeq(dds)
res <- results(dds)

#res <- results(dds, contrast=c("treatment","Phosphate deficiency", "Phosphate sufficiency"))
res <- results(dds, contrast=c("treatment","3hrs red light irradiation (7Ã«_molm-2s-1)", "none"))
resUpregulated <- subset(res, log2FoldChange > 0 & padj < 0.05)
resDownregulated <- subset(res, log2FoldChange < 0 & padj < 0.05)
