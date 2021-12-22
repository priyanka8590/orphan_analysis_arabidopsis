import os, sys

normalized_types_filename = sys.argv[1]
sra_studies_filename = sys.argv[2]
output_filename = sys.argv[3]

normalized_types_file = open(normalized_types_filename, "r")
sra_studies_file = open(sra_studies_filename, "r")
output_file = open(output_filename, "w")

normalized_types_list = []
for line in normalized_types_file:
    normalized_type = list.strip().split("\n")
    normalized_types_list.append(normalized_type)
    
sra_studies_list = []
for line_2 in sra_studies_file:
    sra_studies = line_2.strip().split("\n")
    sra_studies_list.append(sra_studies)
    
for normalized_types in normalized_types_list:
    file_to_be_created = "merged_"+normalized_types+".tsv"
    filename_to_be_created = open(file_to_be_created, "w")
    for sra_study in sra_studies_list:
        filename_to_be_used = sra_study+"_counts_for_norm_within_studies.tsv"
        sra_study_gene_to_counts = {}
        file_to_be_used = open(filename_to_be_used, "r")
        for line in file_to_be_used:
            gene = line.strip().split("\t")[0]
            counts = line.strip().split("\t")[1:]
            if gene not in sra_study_gene_to_counts:
                sra_study_gene_to_counts[gene] = counts
    for sra_study_in_dict in sra_study_gene_to_counts:
        filename_to_be_created.write(sra_study_gene_to_counts+"\t"+"\t")
        