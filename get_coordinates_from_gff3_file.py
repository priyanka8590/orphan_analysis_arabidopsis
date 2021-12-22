import os, sys

gff3_filename = sys.argv[1]
gene_list_filename = sys.argv[2]
output_filename = sys.argv[3]

gff3_file = open(gff3_filename, "r")
gene_list_file = open(gene_list_filename, "r")
output_file = open(output_filename, "w")

id_to_coordinates = {}
for line in gff3_file:
    if line[0] == "#":continue
    chromosome, source, typeofgene, start, end, score, strand, phase, attributes = line.strip().split("\t")
    if typeofgene != "CDS":continue
    parent_id = attributes.strip().split(";")[1].split("Parent=")[-1]
    if parent_id not in id_to_coordinates:
        id_to_coordinates[parent_id] = [chromosome, start, end]
        
gff3_file.close()
        
ids_from_gene_list = []
first_line = gene_list_file.readline()
for line_2 in gene_list_file:
    qseqid = line_2.strip().split()[1][1:-1]
    ids_from_gene_list.append(qseqid)
    
#print(ids_from_gene_list)

gene_list_file.close()


for gene_id in id_to_coordinates:
    if gene_id in ids_from_gene_list:
        #print (gene_id)
        #print(id_to_coordinates[gene_id][0], id_to_coordinates[gene_id][1], id_to_coordinates[gene_id][2])
        output_file.write(id_to_coordinates[gene_id][0]+"\t"+id_to_coordinates[gene_id][1]+"\t"+id_to_coordinates[gene_id][2]+"\t"+gene_id+":"+id_to_coordinates[gene_id][0]+":"+id_to_coordinates[gene_id][1]+"-"+id_to_coordinates[gene_id][2])
        output_file.write("\n")
        