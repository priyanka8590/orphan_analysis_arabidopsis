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

setwd("~/Box/analysis_orphan_gene_aim_3/deseq_analysis/cold_stress/Exp1_dex_col0")
setwd("~/Box/analysis_orphan_gene_aim_3")
salmon_counts <- read.table(file = "salmon_consolidated_all_5210.tsv", sep = "\t", header = TRUE)
count_matrix <- read.table("final_expression_table_with_high_depth_25.tsv", sep = "\t", header = TRUE)
studies_with_one_run <- read.delim("runs_to_be_removed")
studies_with_one_run <- as.data.frame(studies_with_one_run)
new_variable <- as.vector(studies_with_one_run$runs)
df <- count_matrix[,!(names(count_matrix) %in% new_variable)]
metadata_table <- read.delim("final_metadata_for_4741_runs_for_deseq2_original.txt")
metadata_table <- metadata_table[!(metadata_table$run_accession %in% new_variable),]
metadata_table_og <- read.delim("final_metadata_for_5210_runs")

salmon_counts <- subset(salmon_counts, Gene %in% df$rowname) #to get filtered genes - 53052 genes , df is final_expression file dataframe
salmon_counts <- as.data.frame(salmon_counts)
salmon_counts_pi_def <- salmon_counts %>% select(Gene, SRR6048692, SRR6048693, SRR6048694, SRR6048695, SRR6048695, SRR6048696, SRR6048697, SRR6048698, SRR6048699, SRR6048700, SRR6048701, SRR6048702, SRR6048703, SRR6048704, SRR6048705, SRR6048706, SRR6048707)
salmon_counts_pi_def <- salmon_counts %>% select(Gene, starts_with('SRR60292'))
salmon_counts_pi_def_mat <- data.matrix(salmon_counts_pi_def)
rownames(salmon_counts_pi_def_mat) <- salmon_counts_pi_def[,1]
salmon_counts_pi_def_mat <- salmon_counts_pi_def_mat[,-1]
salmon_counts_pi_def_mat <- round(salmon_counts_pi_def_mat)

#To get metadata
library(dplyr)
#dea_anno <- filter(metadata_table, study_accession %in% 'SRP068145')
dea_exp1_anno <- filter(metadata_table, study_accession %in% 'SRP117933')
dea_exp1_anno <- read.delim("metadata_table_cold_stress_1")
dea_exp1_anno <- as.data.frame(dea_exp1_anno)
dea_exp1_anno<-mutate(dea_exp1_anno,time=as.character(time))
dea_exp1_anno<-mutate(dea_exp1_anno,temperature=as.character(temperature))
dea_exp1_anno<-mutate(dea_exp1_anno,time=sapply(strsplit(dea_exp1_anno$time, split='time - ', fixed=TRUE),function(x) (x[2])))
dea_exp1_anno<-mutate(dea_exp1_anno,temperature=sapply(strsplit(dea_exp1_anno$temperature, split='temp - ', fixed=TRUE),function(x) (x[2])))
dea_exp1_anno$treatment <- c("22c", "27c", "22c", "22c", "27c", "22c", "27c", "22c", "27c", "22c", "27c", "22c", "27c", "22c", "27c", "27c")
dea_exp1_anno$time <- c("22", "16", "4", "20", "4", "0", "16", "1", "12", "8", "22", "8", "1", "0", "20", "12")
rownames(dea_exp1_anno) <- dea_exp1_anno[,1]
dea_exp1_anno <- dea_exp1_anno[,-1]
dea_exp1_anno <- dea_exp1_anno[order(row.names(dea_exp1_anno)),]
dea_exp1_anno <- dea_exp1_anno[,c("time", "treatment")]
dea_exp1_anno$treatment <- factor(dea_exp1_anno$treatment)
dea_exp1_anno$time <- factor(dea_exp1_anno$time)

col_order <- c("Gene", rownames(dea_exp1_anno))
salmon_counts_pi_def <- salmon_counts_pi_def[,col_order]

dds <- DESeqDataSetFromMatrix(countData = salmon_counts_pi_def_mat, colData = dea_exp1_anno, design = ~ treatment + time + treatment:time)

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
res <- results(dds, contrast = c("treatment", "25c", "16c"))
resUpregulated <- subset(res, log2FoldChange > 0 & padj < 0.05)
resDownregulated <- subset(res, log2FoldChange < 0 & padj < 0.05)
write.table(resUpregulated, file = "upregulated_in treatment.tsv", sep = "\t", quote = F, row.names = TRUE)
write.table(resDownregulated, file = "downregulated_in treatment.tsv", sep = "\t", quote = F, row.names = TRUE)

res <- results(dds, contrast = list(c("treatment_25c_vs_16c", "treatment25c.genotypesdg261")))
resUpregulated <- subset(res, log2FoldChange > 0 & padj < 0.05)
resDownregulated <- subset(res, log2FoldChange < 0 & padj < 0.05)
write.table(resUpregulated, file = "upregulated_in_treatment_for_genotype_3.tsv", sep = "\t", quote = F, row.names = TRUE)
write.table(resDownregulated, file = "downregulated_in_treatment_for_genotype_3.tsv", sep = "\t", quote = F, row.names = TRUE)

res <- results(dds, contrast = list(c("treatment_25c_vs_16c", "treatment25c.genotypesdg82")))
resUpregulated <- subset(res, log2FoldChange > 0 & padj < 0.05)
resDownregulated <- subset(res, log2FoldChange < 0 & padj < 0.05)
write.table(resUpregulated, file = "upregulated_in_treatment_for_genotype_2.tsv", sep = "\t", quote = F, row.names = TRUE)
write.table(resDownregulated, file = "downregulated_in_treatment_for_genotype_2.tsv", sep = "\t", quote = F, row.names = TRUE)

res <- results(dds, contrast = c(0, -1, 1, 0, 0, 0))
resUpregulated <- subset(res, log2FoldChange > 0 & padj < 0.05)
resDownregulated <- subset(res, log2FoldChange < 0 & padj < 0.05)
write.table(resUpregulated, file = "upregulated_in_treatment_for_genotype_2_vs_genotype_3.tsv", sep = "\t", quote = F, row.names = TRUE)
write.table(resDownregulated, file = "downregulated_in_treatment_for_genotype_2_vs_genotype_3.tsv", sep = "\t", quote = F, row.names = TRUE)

res <- results(dds, name = "treatment25c.genotypesdg261")
resUpregulated <- subset(res, log2FoldChange > 0 & padj < 0.05)
resDownregulated <- subset(res, log2FoldChange < 0 & padj < 0.05)
write.table(resUpregulated, file = "upregulated_in_interaction_for_genotype_3_vs_genotype_1.tsv", sep = "\t", quote = F, row.names = TRUE)
write.table(resDownregulated, file = "downregulated_in_interaction_for_genotype_3_vs_genotype_1.tsv", sep = "\t", quote = F, row.names = TRUE)

res = results(dds, contrast = list("treatment25c.genotypesdg261", "treatment25c.genotypesdg82"))
resUpregulated <- subset(res, log2FoldChange > 0 & padj < 0.05)
resDownregulated <- subset(res, log2FoldChange < 0 & padj < 0.05)
write.table(resUpregulated, file = "upregulated_in_interaction_for_genotype_3_vs_genotype_2.tsv", sep = "\t", quote = F, row.names = TRUE)
write.table(resDownregulated, file = "downregulated_in_interaction_for_genotype_3_vs_genotype_2.tsv", sep = "\t", quote = F, row.names = TRUE)