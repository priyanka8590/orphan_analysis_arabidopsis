import os, sys

normalized_expression_filename = sys.argv[1]
gene_metadata_filename = sys.argv[2]
final_filename_for_mog = sys.argv[3]

normalized_expression_file = open(normalized_expression_filename, 'r')
gene_metadata_file = open(gene_metadata_filename, 'r')
final_file_for_mog = open(final_filename_for_mog, 'w')

run_to_counts = {}
runs = []
first_line = normalized_expression_file.readline()
for line in normalized_expression_file:
    run = line.strip().split("\t")[0]
    counts = line.strip().split("\t")[1:]
    runs.append(run)
    run_to_counts[run] = counts
    
run_to_metadata = {}
first_line_metadata_file = gene_metadata_file.readline()
for line_2 in gene_metadata_file:
    run = line_2.strip().split("\t")[0]
    print(run)
    run = run.replace('"', '')
    metadata = line_2.strip().split("\t")[1:]
    run_to_metadata[run] = metadata
    
#final_file_for_mog.write("ID"+"\t"+first_line_metadata_file+"\t"+first_line)
#final_file_for_mog.write("\n")
for run in run_to_counts:
    if run in runs:
        final_file_for_mog.write(run+"\t"+"\t".join(run_to_metadata[run])+"\t"+"\t".join(run_to_counts[run]))
        final_file_for_mog.write("\n")
        
normalized_expression_file.close()
gene_metadata_file.close()
final_file_for_mog.close()
        
    
    