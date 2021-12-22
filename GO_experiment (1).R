library(plyr)
library(clusterProfiler)
#This script is for writing GO analysis result of each cluster for experiment data
newcluster <- read.csv("JL_2018-11-12_scnorm25cluster.csv",header = T)
newcluster <- newcluster[grepl("Y|Q",newcluster$qseqid),]
count <- as.data.frame(table(newcluster[,1]))
cluster6<-as.character(count$Var1[count$Freq>=6])
all <- read.csv("all.csv",header = T)$gene
GO <- read.csv("GO_annotation.csv", header = F)
BP <- subset(GO, V5=="biological_process")
MF <- subset(GO, V5=="molecular_function")
CC <- subset(GO, V5=="cellular_component")
term2gene_MF <- data.frame(MF$V3,MF$V2)
term2name_MF <- data.frame(MF$V3,MF$V4)
term2gene_CC <- data.frame(CC$V3,CC$V2)
term2name_CC <- data.frame(CC$V3,CC$V4)
term2gene_BP <- data.frame(BP$V3,BP$V2)
term2name_BP <- data.frame(BP$V3,BP$V4)

for (i in cluster6) {
  cluster_genes <- as.character(newcluster$qseqid[newcluster[,1]==i])
  
  new_BP <- enricher(gene = cluster_genes,
                     minGSSize = 1,
                     pvalueCutoff  = 1, 
                     pAdjustMethod = "BH", 
                     universe = all, 
                     TERM2GENE = term2gene_BP, 
                     TERM2NAME = term2name_BP)
  if(!is.null(new_BP)){
    result_BP <- new_BP@result
    result_BP <- subset(result_BP,p.adjust<0.05)
    if(nrow(result_BP)>0){
      result_BP$ontology <- "BP"
      result_BP$cluster <- i} 
    else {
      result_BP$ontology <- character()
      result_BP$cluster <- integer()}} else {
        result_BP <- data.frame(ID=character(),
                                Description=character(),
                                GeneRatio=character(),
                                BgRatio=character(),
                                pvalue=numeric(),
                                p.adjust=numeric(),
                                qvalue=numeric(),
                                geneID=character(),
                                Count=integer(),
                                ontology=character(),
                                cluster=integer()
        )
      }
  
  new_CC <- enricher(gene = cluster_genes,
                     minGSSize = 1,
                     pvalueCutoff  = 0.05, 
                     pAdjustMethod = "BH", 
                     universe = all, 
                     TERM2GENE = term2gene_CC, 
                     TERM2NAME = term2name_CC)
  if(!is.null(new_CC)){
    result_CC <- new_CC@result
    result_CC <- subset(result_CC,p.adjust<0.05)
    if(nrow(result_CC)>0){
      result_CC$ontology <- "CC"
      result_CC$cluster <- i} 
    else {
      result_CC$ontology <- character()
      result_CC$cluster <- integer()}} else {
        result_CC <- data.frame(ID=character(),
                                Description=character(),
                                GeneRatio=character(),
                                BgRatio=character(),
                                pvalue=numeric(),
                                p.adjust=numeric(),
                                qvalue=numeric(),
                                geneID=character(),
                                Count=integer(),
                                ontology=character(),
                                cluster=integer()
        )
      }
  
  new_MF <- enricher(gene = cluster_genes,
                     minGSSize = 1,
                     pvalueCutoff  = 0.05, 
                     pAdjustMethod = "BH", 
                     universe = all, 
                     TERM2GENE = term2gene_MF, 
                     TERM2NAME = term2name_MF)
  if(!is.null(new_MF)){
    result_MF <- new_MF@result
    result_MF <- subset(result_MF,p.adjust<0.05)
    if(nrow(result_MF)>0){
      result_MF$ontology <- "MF"
      result_MF$cluster <- i} 
    else {
      result_MF$ontology <- character()
      result_MF$cluster <- integer()}} else {
        result_MF <- data.frame(ID=character(),
                                Description=character(),
                                GeneRatio=character(),
                                BgRatio=character(),
                                pvalue=numeric(),
                                p.adjust=numeric(),
                                qvalue=numeric(),
                                geneID=character(),
                                Count=integer(),
                                ontology=character(),
                                cluster=integer()
        )
      }
  
  result_i <- rbind(result_BP,result_CC,result_MF)
  write.table(
    result_i,
    paste("out_newcombined/",i,".txt",sep = ""),
    sep       = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote     = FALSE)
}