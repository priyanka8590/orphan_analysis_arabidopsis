import os, sys

all_sample_metadata_filename=sys.argv[1]
all_samples_filename=sys.argv[2]
output_filename=sys.argv[3]

all_samples_metadata_file=open(all_sample_metadata_filename,'r',encoding="ISO-8859-1")
all_samples_file=open(all_samples_filename,'r')
output_file=open(output_filename,'w')

run_accession_to_metadata={}
first_line=all_samples_metadata_file.readline()
for line in all_samples_metadata_file:
    run_accession=line.strip().split("\t")[0]
    all_metadata=line.strip().split("\t")[1:]
    if run_accession not in run_accession_to_metadata:
        run_accession_to_metadata[run_accession]=all_metadata
print (run_accession_to_metadata)

all_samples_metadata_file.close()
samples_list=[]
for line in all_samples_file:
    sample=line.strip()
    samples_list.append(sample)
all_samples_file.close()
    
output_file.write(first_line)
"""for run_accession in run_accession_to_metadata:
    if run_accession in samples_list:
        output_file.write(run_accession+"\t"+"\t".join(map(str,run_accession_to_metadata[run_accession])))
        output_file.write("\n")"""
        
for sample in samples_list:
    output_file.write(sample+"\t"+"\t".join(map(str,run_accession_to_metadata[sample])))
    output_file.write("\n")
        
    
    