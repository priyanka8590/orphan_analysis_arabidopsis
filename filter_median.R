library(dplyr)
setwd("/Users/bhandary/Box/analysis_orphan_gene_aim_3/")
setwd("/work/LAS/mash-lab/bhandary/analysis_orphan_gene_aim_3/")
setwd("~/Box/analysis_orphan_gene_aim_3")
expression_file <- read.table(file = "salmon_consolidated_all_5210.tsv", sep = "\t", header = TRUE)
#expression_file <- read.table(file = "final_high_depth_counts.tsv", sep = "\t", header = TRUE)
final_high_depth_counts <- read.csv("~/Box/analysis_orphan_gene_aim_3/final_high_depth_counts.tsv", sep="")
expression_df <- as.data.frame(final_high_depth_counts)
final_df<- expression_df[,-1]
final_df<-as.data.frame(final_df)
rownames(final_df) <- expression_df[,1]
#final_df <- lapply(expression_df, as.numeric)
final_df$Rmean <- rowMeans(final_df)
#final_df <- final_df[final_df$Rmean >= 1,]
#final_df <- filter(final_df, Rmean >= 1)
final_df <- final_df[rowSums(final_df > 0.75 ) >= 25, ]
#final_df[apply(final_df>1,100,any),]
#final_df[apply(final_df, MARGIN = 1, function(x) any(x>1))]
final_df <- tibble::rownames_to_column(final_df)
#rownames(final_df) <- rownames(expression_df)
write.table(final_df, file = "~/Box/analysis_orphan_gene_aim_3/final_expression_table_with_high_depth_25.tsv", sep = "\t", quote = F, row.names = F)


