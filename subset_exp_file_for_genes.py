import os, sys

expression_matrix_filename = sys.argv[1]
gene_filename = sys.argv[2]
subsetted_matrix_filename = sys.argv[3]

expression_matrix_file = open(expression_matrix_filename, "r")
gene_file = open(gene_filename, "r")
subsetted_matrix_file = open(subsetted_matrix_filename, "w")

gene_to_counts = {}
header = expression_matrix_file.readline()
for line in expression_matrix_file:
    gene = line.strip().split("\t")[0]
    counts = line.strip().split("\t")[1:]
    if gene not in gene_to_counts:
        gene_to_counts[gene] = counts

expression_matrix_file.close()
    
ids_from_gene_list = []
first_line = gene_file.readline()
for line_2 in gene_file:
    qseqid = line_2.strip().split()[1][1:-1]
    ids_from_gene_list.append(qseqid)    
    
gene_file.close()

subsetted_matrix_file.write(header)
for gene in gene_to_counts:
    if gene in ids_from_gene_list:
        subsetted_matrix_file.write(gene+"\t"+"\t".join(map(str,gene_to_counts[gene])))
        subsetted_matrix_file.write("\n")
    
