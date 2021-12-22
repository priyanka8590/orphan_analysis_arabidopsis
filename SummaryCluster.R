library(plyr)
library(ggplot2)
library(clusterProfiler)
library(dplyr)
library(data.table)
all <- read.csv("all.csv",header = T)$gene
GO <- read.csv("GO_annotation.csv", header = F)
term2gene <- data.frame(GO$V3,GO$V2)
term2name <- data.frame(GO$V3,GO$V4)
BP <- subset(GO, V5=="biological_process")
MF <- subset(GO, V5=="molecular_function")
CC <- subset(GO, V5=="cellular_component")
term2gene_MF <- data.frame(MF$V3,MF$V2)
term2name_MF <- data.frame(MF$V3,MF$V4)
term2gene_CC <- data.frame(CC$V3,CC$V2)
term2name_CC <- data.frame(CC$V3,CC$V4)
term2gene_BP <- data.frame(BP$V3,BP$V2)
term2name_BP <- data.frame(BP$V3,BP$V4)
file <- c("JL_2018-07-19_top25_cluster.csv","cpm0.6.csv","JL_2018-11-12_scnormannoc.csv","JL_2018-11-12_o25cluster.csv")

for (m in file) {
  positive <- read.csv(m,header = T)
  count <- as.data.frame(table(positive$cluster))
  c6<-subset(count,Freq>=6)
  cluster6<-as.character(c6$Var1)
  set.seed(100)
  sample <- as.data.frame(replicate(100,sample(positive$gene)))
  sample$cluster <- positive$cluster
  bestp_BP_ex <- data.frame(adj.p = numeric(), ontology = character())
  bestp_CC_ex <- data.frame(adj.p = numeric(), ontology = character())
  bestp_MF_ex <- data.frame(adj.p = numeric(), ontology = character())
  for (i in cluster6) {
    cluster_genes <- subset(positive,cluster==i)
    cluster_genes <- as.character(cluster_genes$gene)
    new_BP <- enricher(gene = cluster_genes,
                       minGSSize = 1,
                       pvalueCutoff  = 1, 
                       pAdjustMethod = "BH", 
                       universe = all, 
                       TERM2GENE = term2gene_BP, 
                       TERM2NAME = term2name_BP)
    if(!is.null(new_BP)){
      result_BP <- subset(new_BP@result,Count > 1)
      bestp_BP_ex[i,1]<- min(result_BP$p.adjust)
      bestp_BP_ex$ontology<-"BP"}
    
    new_CC <- enricher(gene = cluster_genes,
                       minGSSize = 1,
                       pvalueCutoff  = 1, 
                       pAdjustMethod = "BH", 
                       universe = all, 
                       TERM2GENE = term2gene_CC, 
                       TERM2NAME = term2name_CC)
    if(!is.null(new_CC)){
      result_CC <- subset(new_CC@result,Count > 1)
      bestp_CC_ex[i,1]<- min(result_CC$p.adjust)
      bestp_CC_ex$ontology<-"CC"}
    
    new_MF <- enricher(gene = cluster_genes,
                       minGSSize = 1,
                       pvalueCutoff  = 1, 
                       pAdjustMethod = "BH", 
                       universe = all, 
                       TERM2GENE = term2gene_MF, 
                       TERM2NAME = term2name_MF)
    if(!is.null(new_MF)){
      result_MF <- subset(new_MF@result,Count > 1)
      bestp_MF_ex[i,1]<- min(result_MF$p.adjust)
      bestp_MF_ex$ontology<-"MF"}
  }
  bestp_BP_ex[bestp_BP_ex == "Inf"] <- NA
  bestp_CC_ex[bestp_CC_ex == "Inf"] <- NA
  bestp_MF_ex[bestp_MF_ex == "Inf"] <- NA
  experiment <- rbind(bestp_BP_ex,bestp_CC_ex,bestp_MF_ex)
  experiment$file <- m
  assign(paste0("experiment_",m),experiment)
  best_p_BP <- data.frame()
  best_p_CC <- data.frame()
  best_p_MF <- data.frame()
  V <- paste0("V",1:100)
  for (j in V) {
    for (i in cluster6) {
      cluster_genes <- subset(sample,cluster==i,select=j)
      colnames(cluster_genes) <- "gene"
      cluster_genes <- as.character(cluster_genes$gene)
      new_BP <- enricher(gene = cluster_genes,
                         minGSSize = 1,
                         pvalueCutoff  = 1, 
                         pAdjustMethod = "BH", 
                         universe = all, 
                         TERM2GENE = term2gene_BP, 
                         TERM2NAME = term2name_BP)
      if(!is.null(new_BP)){
        result_BP <- subset(new_BP@result,Count > 1)
        best_p_BP[i,j]<- min(result_BP$p.adjust)}
      
      new_CC <- enricher(gene = cluster_genes,
                         minGSSize = 1,
                         pvalueCutoff  = 1, 
                         pAdjustMethod = "BH", 
                         universe = all, 
                         TERM2GENE = term2gene_CC, 
                         TERM2NAME = term2name_CC)
      if(!is.null(new_CC)){
        result_CC <- subset(new_CC@result,Count > 1)
        best_p_CC[i,j]<- min(result_CC$p.adjust)}
      
      new_MF <- enricher(gene = cluster_genes,
                         minGSSize = 1,
                         pvalueCutoff  = 1, 
                         pAdjustMethod = "BH", 
                         universe = all, 
                         TERM2GENE = term2gene_MF, 
                         TERM2NAME = term2name_MF)
      if(!is.null(new_MF)){
        result_MF <- subset(new_MF@result,Count > 1)
        best_p_MF[i,j]<- min(result_MF$p.adjust)}
    }
    paste0("finish file ",m," for random set ",j)
  }
  best_p_BP[best_p_BP == "Inf"] <- NA
  best_p_CC[best_p_CC == "Inf"] <- NA
  best_p_MF[best_p_MF == "Inf"] <- NA
  mean_BP<- as.data.frame(apply(best_p_BP,2,function(x)mean(x,na.rm=T)))
  colnames(mean_BP) <- "BP"
  mean_CC <- as.data.frame(apply(best_p_CC,2,function(x)mean(x,na.rm=T)))
  colnames(mean_CC) <- "CC"
  mean_MF <- as.data.frame(apply(best_p_MF,2,function(x)mean(x,na.rm=T)))
  colnames(mean_MF) <- "MF"
  random <- cbind(mean_BP,mean_CC,mean_MF)
  random$file <- m
  assign(paste0("random_",m),random)
}
experiment <- do.call(rbind, lapply( paste0("experiment_", file),get))
experiment_mean <- summarise(group_by(experiment,ontology, file),mean(adj.p,na.rm = T))
colnames(experiment_mean)[3] <- "mean" 
random <- do.call(rbind, lapply( paste0("random_", file),get))
random_t <- melt(random,value.name = "adj.p",variable.name = "ontology")
all <- merge(random_t,experiment_mean,by=c("file","ontology"))
summary <- as.data.frame(summarise(group_by(all,file,ontology),mean(mean),mean(adj.p,na.rm = T),sd(adj.p,na.rm = T),median(adj.p,na.rm = T)))
colnames(summary)[3:6] <- c("ex_mean","mean","sd","median")
attach(summary)
summary$zscore <- abs((ex_mean-mean)/sd)
summary$distance <- abs(ex_mean-median)
write.csv(summary,"JL_2018-01-16_summary.csv",row.names = F)
orphan <- read.csv("orphan.csv")
summary2 <- matrix(nrow=6,ncol = 6)
for(i in 1){
  positive <- read.csv(file[i],header = T)
  ng <- nrow(positive)
  count <- as.data.frame(table(positive$cluster))
  colnames(count)[1] <- "cluster"
  nc <- nrow(count)
  ngmax <- max(count$Freq)
  c6<-subset(count,Freq>=6)
  nc6 <- nrow(c6)
  orinc <- merge(positive,orphan,by="gene")
  norc <- nrow(orinc)
  or_c6 <- merge(c6,orinc,by="cluster")
  norc6 <- nrow(or_c6)
  summary2[i,] <- c(ng,nc,ngmax,nc6,norc,norc6)
}
summary2 <- as.data.frame(summary2)
row.names(summary2) <- file
colnames(summary2) <- c("gene in clusters","cluster","gene in largest cluster","cluster with gene>=6","orphan in clusters","orphan in cluster with gene>=6")
write.csv(summary2,"JL_2018-01-16_summary2.csv",row.names = F)

ggplot(data=random_t,mapping=aes(x=adj.p)) + 
  geom_histogram(binwidth=.001,colour="black",fill="white") + 
  geom_segment(data=summary,aes(x=ex_mean,y=5,xend=ex_mean,yend=0), arrow = arrow(length = unit(0.2, "cm")), colour="red",inherit.aes=FALSE) +
  facet_grid(file~ontology) +
  labs(x = "best p-value", title = "Random Sample for Each Threshold") +
  theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 20))

cpm0.6 <- subset(random_t,file=="cpm0.6.csv")
lab <- c('BP'="Biologcal Process",
         'CC'="Cellular Component",
         'MF'="Molecular Function")
ggplot(data=cpm0.6,mapping=aes(x=adj.p)) + 
  geom_histogram(binwidth=.001,colour="black",fill=c("white")) + 
  geom_segment(data=summary[1:3,],aes(x=ex_mean,y=3,xend=ex_mean,yend=0), arrow = arrow(length = unit(0.2, "cm")), colour="red",inherit.aes=FALSE,size=2.5) +
  facet_grid(ontology~.,labeller = as_labeller(lab)) +
  labs(x = "Best p-value", y = "Count") +
  theme(axis.title.x = element_text(hjust = 0.5,face = "bold",size = 28),
        axis.title.y = element_text(hjust = 0.5,face = "bold",size = 28),
        strip.text.y = element_text(face = "bold",size = 19.5),
        axis.text = element_text(face = "bold",size = 22)) +
  xlim(0,0.125)

ggplot(data=random_t,mapping=aes(x=adj.p)) + 
  geom_histogram(binwidth=.001,colour="black",fill=c("white")) + 
  geom_segment(data=summary[1:3,],aes(x=ex_mean,y=3,xend=ex_mean,yend=0), arrow = arrow(length = unit(0.05, "npc")), colour="red",inherit.aes=FALSE,size=2) +
  facet_grid(ontology~.) +
  labs(x = "Best p-value", y = "Count") +
  theme(axis.title.x = element_text(hjust = 0.5,face = "bold",size = 28),
        axis.title.y = element_text(hjust = 0.5,face = "bold",size = 28),
        strip.text.y = element_text(face = "bold",size = 19.5),
        axis.text = element_text(face = "bold",size = 22)) +
  xlim(0,0.125)
