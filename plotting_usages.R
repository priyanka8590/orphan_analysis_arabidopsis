library(ggplot2)
setwd("~/Box/analysis_orphan_gene_aim_3/")
aa_usage <- read.table("aminoacidusage", header = TRUE, sep = "\t", quote = "\"")
aa_usage <- as.data.frame(aa_usage)
aa_usage <- aa_usage %>% group_by(Type_of_Amino_acid) %>% mutate(ratio = Frequency/Frequency[PS=="Non_Genic_ORFs"])
p <- ggplot(aa_usage, aes(x = Type_of_Amino_acid, y = ratio, group = PS)) + geom_line(aes(colour = PS))
p

aa_usage <- read.table("aminoacidpropusage", header = TRUE, sep = "\t", quote = "\"")
aa_usage <- as.data.frame(aa_usage)
p <- ggplot(aa_usage, aes(x = Amino_acid, y = Frequency, group = PS)) + geom_line(aes(colour = PS))

codon_usage <- read.table("codonusage", header = TRUE, sep = "\t", quote = "\"")
codon_usage <- as.data.frame(codon_usage)
p <- ggplot(codon_usage, aes(x = Codon, y = Frequency, group = PS)) + geom_line(aes(colour = PS)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
