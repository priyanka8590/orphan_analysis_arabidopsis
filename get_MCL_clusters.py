import argparse
import copy
import logging
import multiprocessing
import os
import sys
import time
import pandas as pd



#from ruffus.proxy_logger import *
#from  scipy.stats import spearmanr, pearsonr, kendalltau
#from sklearn import metrics

#import numpy as np
#import pandas as pd
#import scipy.sparse as sparse
#import pprint

def parseCommandLineArguments():
    parser = argparse.ArgumentParser(prog="generate_MCL_clusters.py",description="")
    
    parser.add_argument("--output_directory","-o",help="Enter the name of the directory where all other operations will be performed",required=True)
    parser.add_argument("--gene_counts","-g",help="Enter the name of the gene counts file",required=False)
    parser.add_argument("--correlation_file","-d",help="Enter the file containing the correlation information for the genes",required=True)
    parser.add_argument("--type_of_normalization","-n",help="Enter the type of normalization used in the file", required=True)
    parser.add_argument("--cpu","-p",help="Enter the number of CPU cores to be used",default=1)
    
    return parser.parse_args()

def main():
    commandLineArg=sys.argv
    if len(commandLineArg)==1:
        print("Please use the --help option to get usage information")
    
    options=parseCommandLineArguments()
    performMCLClustering(options)
    
    
def runMCL(eachinput):
    options,threshold,raw_correlation_file,filename_for_mcl_input,filename_for_mcl_output_prefix=eachinput
    if os.path.exists(filename_for_mcl_input)==False:
        fhr=open(raw_correlation_file,"r")
        first_line=fhr.readline()
        fhw=open(filename_for_mcl_input,"w")
        fhw.write("---8<------8<------8<------8<------8<---\n")
        for line in fhr:
            gene1,gene2,correlation=line.strip().split("\t")
            if float(correlation)>threshold:
                fhw.write("\t".join([gene1,gene2,"1"])+"\n")
            else:
                fhw.write("\t".join([gene1,gene2,"0"])+"\n")
        fhw.write("--->8------>8------>8------>8------>8---\n")
        fhw.close()
        fhr.close()
        
    cmd="mcl "+filename_for_mcl_input+" --abc "
    cmd+=" -te "+str(options.cpu)
    cmd+=" -o "+filename_for_mcl_output_prefix+".cluster"
    
    if os.path.exists(filename_for_mcl_output_prefix+".cluster")==False:
        os.system(cmd)
        

def performMCLClustering(options):
    raw_correlation_file = options.correlation_file
    thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
    norm = options.type_of_normalization
    for threshold in thresholds:
        filename_for_mcl_input=options.output_directory+"/spearman_threshold_"+str(threshold)+"_norm_"+str(norm)+".correlation"
        filename_for_mcl_output_prefix=options.output_directory+"/spearman_threshold_"+str(threshold)+"_norm_"+str(norm)
        runMCL([options, float(threshold),raw_correlation_file,filename_for_mcl_input,filename_for_mcl_output_prefix])
    
if __name__ == "__main__":
    main()