import os, sys
import argparse
from pprint import pprint
from sklearn.metrics.cluster import adjusted_rand_score
#import numpy as np
#import pandas as pd

#################################################################################
# Usage
# python find_best_normalization_using_rand_index.py 
# 
#################################################################################


def parseCommandLineArguments():
    parser = argparse.ArgumentParser(prog="find_best_normalization_using_rand_index.py",
                                     description="""Performs all related operations to find best normalization""")
    # Mandatory arguments
    parser.add_argument("--norm_counts_directory","-n",help="Enter the name of the directory with normalized files",required=True)
    parser.add_argument("--output","-o",help="Enter the name of the output file",required=True)
    parser.add_argument("--regulon","-r",help="Enter the name of the regulon file",required=True)
    parser.add_argument("--genes_only_in_new_clusters","-g",help="Enter the name of the file containing genes only found in the new MCL clusters",required=True)
    return parser.parse_args()

def main():
    commandLineArg=sys.argv
    if len(commandLineArg)==1:
        print("Please use the --help option to get usage information")
    options=parseCommandLineArguments()
    regulon_to_genes,list_of_genes_in_regulons,final_list_of_regulons,list_of_genes_only_in_clusters = readGenesToRegulonFile(options)
    checkClusterWithRandIndex(options, regulon_to_genes,list_of_genes_in_regulons,final_list_of_regulons,list_of_genes_only_in_clusters)
    
def readMCLClusterFile(options, norm_type, threshold):
    cluster_file = options.norm_counts_directory+"/spearman_threshold_"+str(threshold)+"_norm_"+str(norm_type)+".cluster"
    dict_of_clusters = {}
    labels_pred=[]
    fhw = open(cluster_file, "r")
    for i,line in enumerate(fhw):
        cluster = line.strip().split("\t")
        #labels_pred.append(i)
        #print (i,cluster)
        dict_of_clusters[i]=cluster
    #print (list_of_clusters)
    #print (type(list_of_clusters))
    fhw.close()
    return dict_of_clusters
    
def readGenesToRegulonFile(options):
    genes_only_in_new_clusters = options.genes_only_in_new_clusters
    fhw=open(genes_only_in_new_clusters, "r")
    list_of_genes_only_in_clusters=[]
    for line in fhw:
        gene=line.strip()
        list_of_genes_only_in_clusters.append(gene)
    fhw.close()
    #print (list_of_genes_only_in_clusters)
    gene_to_regulon_file = options.regulon
    fhw = open(gene_to_regulon_file, "r")
    regulon_to_genes={}
    list_of_genes_in_regulons=[]
    list_of_regulons=[]
    for line in fhw:
        gene, regulon=line.strip().split("\t")
        if gene not in list_of_genes_only_in_clusters:continue
        list_of_genes_in_regulons.append(gene)
        list_of_regulons.append(regulon)
        if regulon not in regulon_to_genes:
            regulon_to_genes[regulon] = []
        regulon_to_genes[regulon].append(gene)
    #print (regulon_to_genes)
    fhw.close()
    final_list_of_regulons = set(list_of_regulons)
    return regulon_to_genes,list_of_genes_in_regulons,final_list_of_regulons,list_of_genes_only_in_clusters

        
def checkClusterWithRandIndex(options, regulon_to_genes,list_of_genes_in_regulons,final_list_of_regulons,list_of_genes_only_in_clusters):
    norm_types = ["cpmnorm", "QNnorm", "normalizedcountsTMM", "qsmooth", "deseq2", "sizenorm"]
    #norm_types = ["deseq2"]
    #thresholds = [0.99,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]
    thresholds = [0.99,0.9,0.8,0.7,0.6,0.5]
    
    """labels_true = []
    #labels_pred = []
    for i, regulon in enumerate(final_list_of_regulons):
        #print (i, regulon)
        for gene in regulon_to_genes[regulon]:
            #print (gene,i,regulon) #9258 genes
            labels_true.append(i)"""
    #print (len(labels_true))
    for norm_type in norm_types:
        for threshold in thresholds:
            dict_of_clusters = readMCLClusterFile(options, norm_type, threshold)
            #print (dict_of_clusters)
            labels_true = []
            for i, regulon in enumerate(final_list_of_regulons):
                #print (i, regulon)
                for gene in regulon_to_genes[regulon]:
                    for cluster_number in dict_of_clusters:
                        if gene in dict_of_clusters[cluster_number]:
                            #print (gene,i,regulon) #9258 genes
                            labels_true.append(i)
            labels_pred = []
            for gene in list_of_genes_only_in_clusters:
                for cluster_number in dict_of_clusters:
                    
                    if gene in dict_of_clusters[cluster_number]:
                        #print (gene, cluster_number,norm_type,threshold)
                        labels_pred.append(cluster_number)
            #print (norm_type, labels_pred)
            #print (len(labels_pred))   
            ari_score = adjusted_rand_score(labels_true, labels_pred)
            print (norm_type, threshold, ari_score)
            """for cluster_number in dict_of_clusters:
                for gene in dict_of_clusters[cluster_number]:
                    if gene in list_of_genes_in_regulons:
                        
                        labels_pred.append(cluster_number)
                        #print (gene, cluster_number)
                        #ari_score = adjusted_rand_score(labels_true, labels_pred)
                        #print (norm_type, threshold, ari_score)"""
    

    
if __name__ == "__main__":
    main()