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

setwd("~/Box/analysis_orphan_gene_aim_3/deseq_analysis/Heat_stress/Exp2_heat_stress_stch4_col0")
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
#salmon_counts_pi_def <- salmon_counts %>% select(Gene, SRR4010821, SRR4010822, SRR4010823, SRR4010824, SRR4010825, SRR4010826, SRR4010827, SRR4010828, SRR4010829, SRR4010830, SRR4010831, SRR4010832, SRR4010833, SRR4010834, SRR4010835, SRR4010836, SRR4010837, SRR4010838)
salmon_counts_pi_def <- salmon_counts %>% select(Gene, starts_with('SRR105745'))
salmon_counts_pi_def_mat <- data.matrix(salmon_counts_pi_def)
rownames(salmon_counts_pi_def_mat) <- salmon_counts_pi_def[,1]
salmon_counts_pi_def_mat <- salmon_counts_pi_def_mat[,-1]
salmon_counts_pi_def_mat <- round(salmon_counts_pi_def_mat)

#To get metadata
library(dplyr)
#dea_anno <- filter(metadata_table, study_accession %in% 'SRP068145')
#dea_exp1_anno <- filter(metadata_table, study_accession %in% 'SRP082470')
dea_exp1_anno <- filter(metadata_table, study_accession %in% 'SRP234442')
rownames(dea_exp1_anno) <- dea_exp1_anno[,1]
dea_exp1_anno <- dea_exp1_anno[,-1]
dea_exp1_anno <- dea_exp1_anno[order(row.names(dea_exp1_anno)),]
#drop <- c("treatment", "genotype")
#dea_exp1_anno <- dea_exp1_anno[,!(names(dea_exp1_anno) %in% drop)]
dea_exp1_anno$treatment <- c("control", "control", "4 degrees for 4 hours", "4 degrees for 4 hours", "4 degrees for 4 hours", "24 degrees for 4 hours", "24 degrees for 4 hours", "24 degrees for 4 hours", "control", "control", "control", "4 degrees for 4 hours", "4 degrees for 4 hours", "4 degrees for 4 hours", "24 degrees for 4 hours", "24 degrees for 4 hours", "24 degrees for 4 hours")
dea_exp1_anno$genotype <- c("Col0", "Col0", "Col0", "Col0", "Col0", "Col0", "Col0", "Col0", "stch4", "stch4", "stch4", "stch4", "stch4", "stch4", "stch4", "stch4", "stch4")
dea_exp1_anno <- dea_exp1_anno[,c("treatment", "genotype")]
#or alternatively,if metadata table doesn't have some runs

dea_exp1_anno$genotype <- factor(dea_exp1_anno$genotype)
dea_exp1_anno$treatment <- factor(dea_exp1_anno$treatment)
col_order <- c("Gene", rownames(dea_exp1_anno))
salmon_counts_pi_def_mat <- salmon_counts_pi_def_mat[,col_order]
dds <- DESeqDataSetFromMatrix(countData = salmon_counts_pi_def_mat, colData = dea_exp1_anno, design = ~ treatment + genotype + treatment:genotype)

#keep <- rowSums(counts(dds)) >= 10
salmon_counts_pi_def_mat[salmon_counts_pi_def_mat > 0] <- floor(salmon_counts_pi_def_mat[salmon_counts_pi_def_mat > 0]) #first remove genes that have no counts in any samples
#salmon_counts_pi_def_mat <- salmon_counts_pi_def_mat[rowSums(salmon_counts_pi_def_mat) > 1,] #never use this
#dds$treatment <- factor(dds$treatment, levels = c("Phosphate deficiency", "Phosphate sufficiency"))
#dds$treatment <- factor(dds$treatment, levels = c("3hrs red light irradiation (7Ã«_molm-2s-1)", "none"))
#dds$treatment <- factor(dds$treatment, levels = c("Phosphate deficiency", "Phosphate sufficiency"))
dds$treatment.n <- as.factor(rep(rep(1:3, each = 3), 2))
dea_exp1_anno$treatment.n <- as.factor(rep(rep(1:3, each = 3), 2))
ml <- model.matrix(~genotype + genotype:treatment.n + genotype:treatment, dea_exp1_anno)

unname(ml)
all.zero <- apply(ml, 2, function(x) all(x==0))
all.zero
idx <- which(all.zero)
ml <- ml[,-idx]
unname(ml)
dds <- DESeq(dds)
dds <- DESeq(dds, full = design(ml), betaPrior = FALSE)
res <- results(dds)
#res <- results(dds, contrast=c("treatment","Phosphate deficiency", "Phosphate sufficiency"))
#res <- results(dds, contrast=c("treatment","3hrs red light irradiation (7Ã«_molm-2s-1)", "none"))
res <- results(dds, contrast = c("treatment", "control", "4 degrees for 4 hours"))
resUpregulated <- subset(res, log2FoldChange > 0 & padj < 0.05)
resDownregulated <- subset(res, log2FoldChange < 0 & padj < 0.05)
write.table(resUpregulated, file = "upregulated_in treatment.tsv", sep = "\t", quote = F, row.names = TRUE)
write.table(resDownregulated, file = "downregulated_in treatment.tsv", sep = "\t", quote = F, row.names = TRUE)


res <- results(dds, contrast = list(c("treatment_4.degrees.for.4.hours_vs_24.degrees.for.4.hours", "treatment4.degrees.for.4.hours.genotypestch4")))
resUpregulated <- subset(res, log2FoldChange > 0 & padj < 0.05)
resDownregulated <- subset(res, log2FoldChange < 0 & padj < 0.05)
write.table(resUpregulated, file = "upregulated_in_24_degrees_for_genotype_stch4.tsv", sep = "\t", quote = F, row.names = TRUE)
write.table(resDownregulated, file = "downregulated_in_24_degrees_for_genotype_stch4.tsv", sep = "\t", quote = F, row.names = TRUE)

res <- results(dds, contrast = list(c("treatment_control_vs_24.degrees.for.4.hours", "treatmentcontrol.genotypestch4")))
resUpregulated <- subset(res, log2FoldChange > 0 & padj < 0.05)
resDownregulated <- subset(res, log2FoldChange < 0 & padj < 0.05)
write.table(resUpregulated, file = "upregulated_in_genotype_control_stch4.tsv", sep = "\t", quote = F, row.names = TRUE)
write.table(resDownregulated, file = "downregulated_in_genotype_control_stch4.tsv", sep = "\t", quote = F, row.names = TRUE)

res <- results(dds, contrast = c("treatment", "control", "24 degrees for 4 hours"))
resUpregulated <- subset(res, log2FoldChange > 0 & padj < 0.05)
resDownregulated <- subset(res, log2FoldChange < 0 & padj < 0.05)
write.table(resUpregulated, file = "upregulated_in_control_vs_24 degrees_for_4 hours.tsv", sep = "\t", quote = F, row.names = TRUE)
write.table(resDownregulated, file = "downregulated_in_control_vs_24 degrees_for_4 hours.tsv", sep = "\t", quote = F, row.names = TRUE)

res <- results(dds, contrast = c(0, -1, 1, 0, 0, 0))
resUpregulated <- subset(res, log2FoldChange > 0 & padj < 0.05)
resDownregulated <- subset(res, log2FoldChange < 0 & padj < 0.05)
write.table(resUpregulated, file = "upregulated_in_treatment_for_genotype_2_vs_genotype_3.tsv", sep = "\t", quote = F, row.names = TRUE)
write.table(resDownregulated, file = "downregulated_in_treatment_for_genotype_2_vs_genotype_3.tsv", sep = "\t", quote = F, row.names = TRUE)

res <- results(dds, name = "treatment4.degrees.for.4.hours.genotypestch4")
resUpregulated <- subset(res, log2FoldChange > 0 & padj < 0.05)
resDownregulated <- subset(res, log2FoldChange < 0 & padj < 0.05)
write.table(resUpregulated, file = "upregulated_in_interaction.tsv", sep = "\t", quote = F, row.names = TRUE)
write.table(resDownregulated, file = "downregulated_in_interaction.tsv", sep = "\t", quote = F, row.names = TRUE)
