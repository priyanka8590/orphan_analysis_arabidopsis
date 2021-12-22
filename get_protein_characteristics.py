import os, sys
import argparse
from collections import Counter

def makeAListOfProteinSequences(options):
    prot_seq = []
    gene_to_peptide_sequence = {}
    with open(options.protein_file, 'r') as f:
        for line in f:
            if line[0] == '>': 
                eb_id = line.strip().split()[0].strip(">")[0:-2].strip('"')
                gene_to_peptide_sequence[eb_id] = f.readline().strip().strip('*')
                
        #line = line.strip().strip('*')
        #prot_seq.append(line)
    #print (gene_to_peptide_sequence)
    #print (gene_to_peptide_sequence["AT1G64633.1"])
    for gene in gene_to_peptide_sequence:
        prot_seq.append(gene_to_peptide_sequence[gene])
    #print (prot_seq)
    return gene_to_peptide_sequence, prot_seq
    
def getAaUsage(gene_to_peptide_dict, ps_to_gene_dict, options):
    aa_usage_dict = {}
    acidic_aa = ['D', 'E']
    basic_aa = ['R', 'K', 'H']
    polar_aa = ['N', 'S', 'T', 'Y', 'C', 'Q']
    non_polar_aa = ['A', 'V', 'L', 'G', 'I', 'M', 'W', 'F', 'P']
    aa_list = ['D', 'E', 'R', 'K', 'H', 'N', 'S', 'T', 'Y', 'C', 'Q', 'A', 'V', 'L', 'G', 'I', 'M', 'W', 'F', 'P']
    output_file = open(options.output_directory+"/aminoacidusage", 'w')
    outputfile_prop = open(options.output_directory+"/aminoacidpropusage", 'w')
    output_file.write("PS"+"\t"+"Type_of_Amino_acid"+"\t"+"Frequency")
    outputfile_prop.write("PS"+"\t"+"Amino_acid"+"\t"+"Type_of_Amino_acid"+"\t"+"Frequency")
    output_file.write("\n")
    outputfile_prop.write("\n")
    for ps in ps_to_gene_dict:
        total_aa = 0
        acidic_aas = []
        basic_aas = []
        polar_aas = []
        non_polar_aas = []
        for aa in aa_list:
            counter = 0
            for gene in ps_to_gene_dict[ps]:
                peptide_seq = gene_to_peptide_dict[gene]
                #print (aa, peptide_seq, peptide_seq.count(aa))
                if aa in acidic_aa:
                    acidic_aas.append(aa)
                elif aa in basic_aa:
                    basic_aas.append(aa)
                elif aa in polar_aa:
                    polar_aas.append(aa)
                elif aa in non_polar_aa:
                    non_polar_aas.append(aa)
                counter = counter + peptide_seq.count(aa)
            if aa not in aa_usage_dict:
                aa_usage_dict[aa] = counter
            #print (ps, aa, counter)
            
            freqDictacidic = dict(Counter(acidic_aas))
            freqDictbasic = dict(Counter(basic_aas))
            freqDictpolar = dict(Counter(polar_aas))
            freqDictnonpolar = dict(Counter(non_polar_aas))
            total_aa = sum(freqDictacidic.values()) + sum(freqDictbasic.values()) + sum(freqDictpolar.values()) + sum(freqDictnonpolar.values())
            total_aa_acidic = sum(freqDictacidic.values())
            total_aa_basic = sum(freqDictbasic.values())
            total_aa_polar = sum(freqDictpolar.values())
            total_aa_non_polar = sum(freqDictnonpolar.values())
            
        for aa in aa_usage_dict:
            output_file.write(ps+"\t"+aa+"\t"+str(aa_usage_dict[aa]/total_aa)+"\n")
        for aa in freqDictacidic:
            outputfile_prop.write(ps+"\t"+aa+"\t"+"acidic"+"\t"+str(round(freqDictacidic[aa]/total_aa_acidic, 2))+"\n")
        for aa in freqDictbasic:
            outputfile_prop.write(ps+"\t"+aa+"\t"+"basic"+"\t"+str(round(freqDictbasic[aa]/total_aa_basic, 2))+"\n")
        for aa in freqDictpolar:
            outputfile_prop.write(ps+"\t"+aa+"\t""polar"+"\t"+str(round(freqDictpolar[aa]/total_aa_polar, 2))+"\n")
        for aa in freqDictnonpolar:
            outputfile_prop.write(ps+"\t"+aa+"\t"+"non_polar"+"\t"+str(round(freqDictnonpolar[aa]/total_aa_non_polar, 2))+"\n")
    #print (aa_usage_dict)
        
def getCodonUsage(options, ps_to_gene_dict):
    #print (ps_to_gene_dict)
    gene_to_cdna_sequence = {}
    output_file_codon_usage = open(options.output_directory+"/codonusage", 'w')
    output_file_codon_usage.write("PS"+"\t"+"Codon"+"\t"+"Frequency")
    output_file_codon_usage.write("\n")
    with open(options.cdna_file, 'r') as f:
        for line in f:
            if line[0] == '>': 
                id = line.strip().split()[0].strip(">").strip('"')
                gene_to_cdna_sequence[id] = f.readline().strip()
    #print (gene_to_cdna_sequence)
    #codon_list = ['UUU', 'UUA', 'UUC', 'UUG', 'CUU', 'CUA', 'CUC', 'CUG', 'AUU', 'AUA', 'AUC', 'AUG', 'GUU', 'GUA', 'GUC', 'GUG', 'UCU', 'UCC', 'UCA', 'UCG', 'CCU', 'CCA', 'CCC', 'CCG', 'ACU', 'ACA', 'ACC', 'ACG', 'GCU', 'GCA', 'GCC', 'GCG', 'UAU', 'UAA', 'UAC', 'UAG', 'CAU', 'CAA', 'CAC', 'CAG', 'AAU', 'AAA', 'AAC', 'AAG', 'GAU', 'GAA', 'GAC', 'GAG', 'UGU', 'UGA', 'UGC', 'UGG', 'CGU', 'CGA', 'CGC', 'CGG','AGU', 'AGA', 'AGC', 'AGG', 'GGU', 'GGA', 'GGC', 'GGG']
    #print (len(codon_list))
    for ps in ps_to_gene_dict:
        tmp_list = []
        for gene in ps_to_gene_dict[ps]:
            cdna_sequence = gene_to_cdna_sequence[gene]
            for codon in range(0, len(cdna_sequence) - 2, 3):
                tmp_list.append(cdna_sequence[codon:codon+3])
        #print (ps, tmp_list)
        freqDict = dict(Counter(tmp_list))
        totalCodons = sum(freqDict.values())
        #print (ps, freqDict, totalCodons)
        for codons in freqDict:
            #print (ps+"\t"+codons+"\t"+str(round(freqDict[codons]/totalCodons, 2)))
            output_file_codon_usage.write(ps+"\t"+codons+"\t"+str(round(freqDict[codons]/totalCodons, 2))+"\n")
        

def getPSInfoForGenes(options):
    ps_to_gene_dict = {}
    ps_info_file = open(options.ps_info_file)
    first_line = ps_info_file.readline()
    for line in ps_info_file:
        gene = line.strip().split("\t")[0].strip('"')
        ps = line.strip().split("\t")[1]
        if ps not in ps_to_gene_dict:
            if gene[0:2] == "AT" and ps == "Arabidopsis thaliana":
                #print (gene, ps)
                ps = "Annotated_Orphans"
            elif gene[0:2]!= "AT":
                #print (gene, ps)
                ps = "Orphan_EB_ORFs"
            ps_to_gene_dict[ps] = []
        ps_to_gene_dict[ps].append(gene)    
    #print (ps_to_gene_dict.keys())
    return ps_to_gene_dict

def getNonGenicOrfs(options):
    non_genic_orfs = []
    with open(options.orfs_file, 'r') as f:
        for line in f:
            if line[0] == '>': 
                orf_id = line.strip().split()[0].strip(">")
                non_genic_orfs.append(f.readline().strip())
    #print (non_genic_orfs)
    return (non_genic_orfs)
        

def parseCommandLineArguments():
    parser = argparse.ArgumentParser(prog="generate_MCL_clusters.py",description="")
    
    parser.add_argument("--output_directory","-o",help="Enter the name of the directory where all other operations will be performed",required=True)
    parser.add_argument("--protein_file","-p",help="Enter the name of the peptide file",required=False)
    parser.add_argument("--ps_info_file","-ps",help="Enter the file containing the phylostratr information for the genes",required=True)
    parser.add_argument("--cdna_file","-c",help="Enter the cdna file", required=True)
    parser.add_argument("--orfs_file","-orf",help="Enter the orfs file outputted by orfipy", required=True)
    parser.add_argument("--cpu","-u",help="Enter the number of CPU cores to be used",default=1)
    
    return parser.parse_args()

def aaAcidUsage(gene_to_peptide_dict, ps_to_gene_dict, options):
    output_file = open(options.output_directory+"/aminoacidusage", 'w')
    outputfile_prop = open(options.output_directory+"/aminoacidpropusage", 'w')
    output_file.write("PS"+"\t"+"Type_of_Amino_acid"+"\t"+"Frequency")
    outputfile_prop.write("PS"+"\t"+"Amino_acid"+"\t"+"Type_of_Amino_acid"+"\t"+"Frequency")
    output_file.write("\n")
    outputfile_prop.write("\n")
    for ps in ps_to_gene_dict:
        aa_list = []
        for gene in ps_to_gene_dict[ps]:
            peptide_seq = gene_to_peptide_dict[gene]
            for aa in peptide_seq:
                aa_list.append(aa)
            #print (aa_list)
        freqDictaa = dict(Counter(aa_list))
        total_aa = sum(freqDictaa.values())
        #print (ps, freqDictaa, total_aa)

        for aa in freqDictaa:
            output_file.write(ps+"\t"+aa+"\t"+str(round(freqDictaa[aa]/total_aa, 2))+"\n")
    return output_file
        
def getAaUsageInNonGenicOrfs(list_of_non_genic_orfs, output_file):
    aa_list = []
    for seq in list_of_non_genic_orfs:
        #counter = 0
        for aa in seq:
            aa_list.append(aa)
    freqDictaa = dict(Counter(aa_list))
    total_aa = sum(freqDictaa.values())
    for aa in freqDictaa:
        output_file.write("Non_Genic_ORFs"+"\t"+aa+"\t"+str(round(freqDictaa[aa]/total_aa, 2))+"\n")
            


def main():
    commandLineArg=sys.argv
    if len(commandLineArg)==1:
        print("Please use the --help option to get usage information")
    
    options=parseCommandLineArguments()
    gene_to_peptide_sequence, prot_seq = makeAListOfProteinSequences(options)
    ps_to_gene_dict = getPSInfoForGenes(options)
    #getAaUsage(gene_to_peptide_sequence, ps_to_gene_dict, options)
    getCodonUsage(options, ps_to_gene_dict)
    output_file = aaAcidUsage(gene_to_peptide_sequence, ps_to_gene_dict, options)
    non_genic_orfs = getNonGenicOrfs(options)
    new_list = list(set(prot_seq).difference(non_genic_orfs))
    print (len(prot_seq))
    print(len(non_genic_orfs))
    print (len(new_list))
    getAaUsageInNonGenicOrfs(new_list, output_file)
    
    
if __name__ == "__main__":
    main()