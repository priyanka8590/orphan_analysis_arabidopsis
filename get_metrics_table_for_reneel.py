import os, sys
import argparse
import numpy as np
import random
#from macpath import join
from numpy.random.mtrand import shuffle
import numpy as np
#import matplotlib.pyplot as plt
# roc curve and auc score
#from sklearn.datasets import make_classification
#from sklearn.neighbors import KNeighborsClassifier
#from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from collections import OrderedDict
#from matplotlib.pyplot import axis
#from sklearn.metrics import plot_roc_curve
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import f1_score
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import roc_auc_score
from sklearn.metrics.ranking import precision_recall_curve

def parseCommandLineArguments():
    """
    Parses the arguments provided through command line.
    Launch python get_metrics_table_for_reneel.py --help for more details
    """
    parser = argparse.ArgumentParser(prog="get_metrics_table_for_mcl.py",description="print this")
    
    parser.add_argument("--datafile","-d",help="Enter the name of the data file",required=True)#d is the data file, microarray or RNA-seq
    parser.add_argument("--output","-o",help="Please enter the name of the output directory",required=True)
    return parser.parse_args()

def main():
    commandLineArg=sys.argv
    if len(commandLineArg)==1:
        print("Please use the --help option to get usage information")
    options = parseCommandLineArguments()
    pca_comps = []
    #print ("I am here!!")
    for i in range(25, 501, 25):
        pca_comps.append(i)
    mcl_thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
    fhw = open("/project/maizegdb/sagnik/bhandary/analysis_regulon_prediction/predict_regulons/genes_to_regulons.tsv", "r")
    regulon_to_gene = {}
    for line in fhw:
        gene, regulon = line.strip().split("\t")
        if regulon not in regulon_to_gene:
            regulon_to_gene[regulon] = []
        regulon_to_gene[regulon].append(gene)
    fhw.close()
    new_file = "/project/maizegdb/sagnik/bhandary/analysis_regulon_prediction/reneel_table_with_metrics.tsv"
    fhw = open(new_file, "w")
    countArray, random_nonphoto_gene = readDataFile(options.datafile)
    for regulon in regulon_to_gene:
        #print ("I am here!!") #final_degree_spearman_425_for_reneel_threshold_numeric_0.3.txt
        for correlation in ["pearson", "spearman"]:
            for pca_comp in pca_comps:
                for mcl_threshold in mcl_thresholds:
                    #print (regulon)
                    #print (correlation)
                    #sys.stdout.flush()
                    mcl_cluster_file = options.output+"/final_degree_"+correlation+"_"+str(pca_comp)+"_for_reneel_threshold_numeric_"+str(mcl_threshold)+".txt"
                    mcl_array, genearray, random_regulon_genes, random_nonregulon_genes = mclClustersWithPathwayInfo(options, mcl_cluster_file, regulon_to_gene[regulon])
                    #print (mcl_array)
                    #sys.stdout.flush()
                    X_training, Y_training = createTrainingDatasetsForMcl(countArray, mcl_array, genearray)
                    #print (Y_training.count(1), Y_training.count(0))
                    #sys.stdout.flush()
                    X_train_initial, X_holdout, Y_train_initial, Y_holdout = train_test_split(X_training, Y_training, test_size=0.10, random_state=42, stratify=Y_training)
                    X_train, X_test, Y_train, Y_test = train_test_split(X_train_initial, Y_train_initial, test_size=0.15, random_state=42, stratify=Y_train_initial)
                    precision, recall, thresholds = makeThresholdForTesting(X_train,Y_train, countArray, genearray, mcl_array, random_regulon_genes,random_nonregulon_genes, X_test, Y_test)
                    if type(precision) == int:
                        line_to_be_written_to_file = f"{precision}\t{recall}\t{thresholds}\t{regulon}\t{pca_comp}\t{mcl_threshold}\t{correlation}"
                        fhw.write(line_to_be_written_to_file)
                        fhw.write("\n")
                        print (line_to_be_written_to_file)
                        sys.stdout.flush()
                    else:
                        for i in range(len(precision)-1):
                            #print ("I am here!!")
                            line_to_be_written_to_file = f"{precision[i]}\t{recall[i]}\t{thresholds[i]}\t{regulon}\t{pca_comp}\t{mcl_threshold}\t{correlation}"
                            fhw.write(line_to_be_written_to_file)
                            fhw.write("\n")
                            print (line_to_be_written_to_file)
                            sys.stdout.flush()
    fhw.close()
    
    
def mclClustersWithPathwayInfo(options, mcl_cluster_file, genes_in_regulon):
    mcl_array = readMCLClusterFile(mcl_cluster_file)
    #genearray = regulon
    random_photo_gene = []
    random_nonphoto_gene = []
    #print (genes_in_regulon)
    #allclusters = list(mcl_array.keys())
    #print (mcl_array.get(allclusters[0]))
    #mcl_photosynthesis_array = {}
    #print(len(genearray))
    for cluster in mcl_array:
        
        for gene in mcl_array[cluster][0]:
            #print (mcl_array[cluster])
            if gene in genes_in_regulon:
                mcl_array[cluster][1].append(1)
            else:
                mcl_array[cluster][1].append(0)
        #print (cluster,len(mcl_array[cluster][0]),mcl_array[cluster][1].count(1)) 
    #print(sum([len(mcl_array[cluster][0]) for cluster in mcl_array])) 
    #for cluster in mcl_array:
        #print (cluster, mcl_array[cluster][1].count(1))
        #sys.stdout.flush()
    return mcl_array, genes_in_regulon, random_photo_gene, random_nonphoto_gene

def readDataFile(datafile):
    countdata = open(datafile, 'r')
    countArray = {}
    first_line = countdata.readline()
    #listofcounts = []
    for line in countdata:
        line = line.strip()
        #print(line)
        line = line.split()
        if ";" in line[0]:
            #line[0] = line[0].split(";")
            for id in line[0].split(";"):
                countArray[id.upper()] = list(map(float,line[1:]))
        #print(line[0])
        else:
            countArray[line[0].upper()] = list(map(float,line[1:]))
            #print(countArray)
            #listofcounts.append(countArray[line[0]])
    #print (countArray)
    genelist = []
    genelist.extend(countArray.keys())
    random.shuffle(genelist)
    random_nonphoto_gene = []
    #print (random_nonphoto_gene)
    return countArray, random_nonphoto_gene

def readMCLClusterFile(mcl_file):
    mcl_clusters = open(mcl_file, 'r')
    mcl_array = {}
    for i, line in enumerate(mcl_clusters):
        line = line.strip().split(",")
        mcl_array[i] = [line,[]]
    #print (mcl_array)
    return (mcl_array)

def createTrainingDatasetsForMcl(countArray, mcl_array, genearray):
    Y_training = []
    X_training = []
    for gene in countArray:
        if gene in genearray:
            Y_training.append(1)
        else:
            Y_training.append(0)
    for gene in countArray.keys():
        X_training.append(gene)
    #print (X_training, Y_training)
    return X_training, Y_training

def makeThresholdForTesting(X_train,Y_train, countArray, genearray, mcl_array, random_photo_genes, random_nonphoto_gene, X_test, Y_test):
    y_true_temp=[]
    for cluster in mcl_array:
        for gene_num,gene in enumerate(mcl_array[cluster][0]):
            if gene in X_test:
                y_true_temp.append(mcl_array[cluster][1][gene_num])
                """print (mcl_array[cluster][1].count(1))
                print (mcl_array[cluster][1].count(0))"""
    #print(y_true_temp.count(1),y_true_temp.count(0))
    p = 0
    #mcl_array_threshold = {}
    #thresholds = {}
    #ground_truth = {ele:0 for ele in random_nonphoto_gene}
    #ground_truth.update({ele:1 for ele in random_photo_genes})
    #fpr = dict()
    #tpr = dict()
    #roc_auc = dict()
    #y_pred = array()
    #y_true = array()
    
    y_pred = [0 for ele in range(len(Y_test))]
    y_true = [0 for ele in range(len(Y_test))]
    y_pred_binary = [0 for ele in range(len(Y_test))]
        #photo_predict = {ele:0 for ele in random_photo_genes + random_nonphoto_gene} 
        #print(photo_predict)
        
    for cluster in mcl_array:
        check_list =[]
        for num, gene in enumerate(mcl_array[cluster][0]):
            if gene in X_test:
                check_list.append("NA")
            else:
                check_list.append(mcl_array[cluster][1][num])
            
        if "NA" not in check_list:continue
        if len(set(check_list))==1 and "NA" in check_list:continue
            #print(check_list,len(mcl_array[cluster][0]))
        p = check_list.count(1)/(len(check_list)-check_list.count("NA"))
                #y_pred.index(check_list.index(gene)) = 1
        #print ("p is ", p, "# genes in training in cluster in regulon ", check_list.count(1), "# genes in training in cluster", (len(check_list)-check_list.count("NA")),"cluster size",len(check_list))
        for ele_num,ele in enumerate(check_list):
            if ele=="NA":
                #print ("MCL cluster", mcl_array[cluster][1][ele_num])
                gene_of_interest=mcl_array[cluster][0][ele_num]
                y_pred[X_test.index(gene_of_interest)]=p
                y_pred_binary[X_test.index(gene_of_interest)]=(1 if p>0.5 else 0) 
                y_true[X_test.index(gene_of_interest)]=mcl_array[cluster][1][ele_num]
    
    #print(y_true.count(1),y_true.count(0))
    #print(y_pred_binary.count(1), y_pred_binary.count(0))
    sys.stdout.flush()
    if y_pred_binary.count(1) < 2 or y_true.count(1) < 2:
        precision = 0
        recall = 0
        thresholds = 0
        """f1_test_score = 0
        matthews_corr = 0
        ROC_score = 0"""
    else:
        precision, recall, thresholds = precision_recall_curve(y_true, y_pred)
        """precision = precision_score(y_true, y_pred, average = 'binary')
        recall = recall_score(y_true, y_pred, average = 'binary')
        f1_test_score = f1_score(y_true, y_pred, average = 'binary')
        matthews_corr = matthews_corrcoef(y_true, y_pred)
        ROC_score = roc_auc_score(y_true, y_pred)"""
        """if len(list(set(y_true)))!=2:
            print (set(y_true))"""
        """for ele_num,ele in enumerate(check_list):
            if ele=="NA":
                gene_of_interest=mcl_array[cluster][0][ele_num]
                y_true[X_test.index(gene_of_interest)]=mcl_array[cluster][1][ele_num]"""
    #print (set(y_true))
        """for ele in y_pred:
            if ele < 0.5:
                y_pred_binary.append(0)
            else:
                y_pred_binary.append(1)"""
    """precision = precision_score(y_true, y_pred_binary, average = 'binary')
    recall = recall_score(y_true, y_pred_binary, average = 'binary')
    f1_test_score = f1_score(y_true, y_pred_binary, average = 'binary')
    matthews_corr = matthews_corrcoef(y_true, y_pred_binary)
    ROC_score = roc_auc_score(y_true, y_pred_binary)"""
    #print (precision+"\t"+recall+"\t"+f1_test_score+"\t"+matthews_corr+"\t"+ROC_score+"\t"+len(X_train)+"\t"+len(X_test))
    #print(precision,Y_train.count(1),Y_test.count(1))
    
    #print("\t".join(list(map(str,[precision,recall,f1_test_score,matthews_corr,ROC_score,len(X_train),len(X_test), Y_train.count(1), Y_test.count(1)]))))
    #print ("\n")
    return precision, recall, thresholds

if __name__ == "__main__":
    main()

