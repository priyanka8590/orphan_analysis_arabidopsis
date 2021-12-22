import os,sys

def mergeDataFileAndGeneMetadata(data_file, gene_metadata, output_file):
    gene_to_counts = {}
    gene_to_metadata = {}
    first_line_data = data_file.readline()
    for line in data_file:
        gene = line.strip().split("\t")[0]
        counts = line.strip().split("\t")[1:]
        gene_to_counts[gene] = counts
    
    first_line_metadata = gene_metadata.readline()
    for line_2 in gene_metadata:
        gene_m = line_2.strip().split("\t")[0]
        metadata = line_2.strip().split("\t")[1:]
        gene_metadata[gene_m] = metadata
    
    output_file.write(gene+"\t"+"\t".join(first_line_metadata)+"\t"+"\t".join(first_line_data))
    output_file.write("\n")
    for gene in gene_to_counts:
        output_file.write(gene+"\t"+"\t".join(gene_metadata[gene])+"\t"+"\t".join(gene_to_counts[gene]))
        output_file.write("\n")

def main():
    data_filename = sys.argv[1]
    gene_metadata_filename = sys.argv[2]
    output_filename = sys.argv[3]
    
    data_file = open(data_filename, 'r')
    gene_metadata_file = open(gene_metadata_filename, 'r')
    output_file = open(output_filename, 'w')
    mergeDataFileAndGeneMetadata(data_file, gene_metadata_file, output_file)
    
    
if __name__ == "__main__":
    main()