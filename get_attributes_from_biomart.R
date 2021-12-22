library(biomaRt)
setwd("~/Box/analysis_orphan_gene_aim_3")
#genes_for_biomart<-read.delim("eb_and_non_db_genes_final.txt")
genes_for_biomart<-read.delim("all_non_eb_genes_25.tsv")
ensembl <- useEnsemblGenomes(biomart = "plants_mart")
datasets <- listDatasets(ensembl)
searchDatasets(mart = ensembl, pattern = "athaliana_eg_gene")
ensembl <- useDataset(dataset = "athaliana_eg_gene", mart = ensembl)

#biomart_obj <- useEnsemblGenomes(biomart = "plants_mart", dataset = "athaliana_eg_gene")
#ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
#                      dataset = "athaliana_eg_gene")
#mart <- useEnsembl(biomart = "ensembl",
#                   dataset = "athaliana_eg_gene")
#ensembl <- useEnsembl(biomart = "genes", dataset = "athaliana_eg_gene")
res<-getBM(attributes = c("external_gene_name","ensembl_gene_id","description","gene_biotype","peptide","chromosome_name", "cds_start","cds_end","ensembl_transcript_id"), filters = "ensembl_transcript_id", values = genes_for_biomart$rowname, mart = ensembl,uniqueRows = TRUE)
res <- as.data.frame(res)
#res<-getBM(attributes = c("ensembl_gene_id","description","gene_biotype","peptide","cds_start","cds_end", "ensembl_transcript_id",), filters = "ensembl_transcript_id", values = genes_for_biomart$rowname, mart = ensembl,uniqueRows = TRUE)
attribute =listAttributes(ensembl)
#write.table(res,file = "metadata_from_biomart_non_eb.txt",sep = "\t", row.names = F)
write.table(res,file = "metadata_from_biomart_all_non_eb_25.tsv",sep = "\t", row.names = F)

#get Jing's delimited file 
genes_from_jing<-read.delim("arabidopsis_all.txt")
write.table(genes_from_jing,file = "arabidopsis_all.tsv",sep = "\t", quote = F, row.names = F)


complete_metadata <- read.delim("complete_gene_metadata_for_non_eb_genes")
colnames(complete_metadata) <- c("transcript_id", "gene_biotype", "peptide_seq", "gene_name", "gene_id", "description", "chromosome", "cds_start", "cds_end", "strata", "chromosome", "transcript_start", "transcript_end", "strand", "cds_len", "cds_gc", "exon_num", "mind", "bind", "dirinf", "maker", "braker", "araport11")
write.table(complete_metadata,file = "complete_metadata_non_eb.txt",sep = "\t", row.names = F)


complete_metadata <- read.delim("complete_metadata_non_eb_and_eb_genes.txt")
complete_metadata <- complete_metadata[,c(4, 6, 10, 2, 17, 16, 15, 3, 7, 12, 13, 8, 9, 14, 1, 18:23)]
write.table(complete_metadata,file = "complete_metadata_non_eb_and_eb_genes.txt",sep = "\t", row.names = F)

tair_description = read.delim("TAIR10_functional_descriptions.txt")
tair_description <- as.data.frame(tair_description)
write.table(tair_description,file = "TAIR_functional_descriptions_tsv.txt",sep = "\t", row.names = F)

run_metadata <- read.delim("final_metadata_for_4741_runs_for_deseq2_original.txt")
write.table(run_metadata,file = "final_4741_metadata.txt",sep = "\t", row.names = F)
