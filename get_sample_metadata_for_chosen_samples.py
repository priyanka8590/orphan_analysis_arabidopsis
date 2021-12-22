import os, sys

#samples_list_filename=sys.argv[1]
smaller_metadata_filename=sys.argv[1]
samples_metadata_filename=sys.argv[2]
final_metadata_filename=sys.argv[3]

#samples_list_file=open(samples_list_filename, 'r')
smaller_metadata_file=open(smaller_metadata_filename,'r')
samples_metadata_file=open(samples_metadata_filename,'r')
final_metadata_file=open(final_metadata_filename,'w')

"""samples_list=[]

for line in samples_list_file:
    sample = list.strip()
    samples_list.append(sample)

samples_list_file.close()"""

samples_list=[]

for line in smaller_metadata_file:
    sample=line.strip().split("\t")[0]
    samples_list.append(sample)

first_line=samples_metadata_file.readline()
sample_to_metadata={}

for line_2 in samples_metadata_file:
    sample=line_2.strip().split("\t")[0]
    metadata=line_2.strip().split("\t")[1:]
    if sample not in sample_to_metadata:
        sample_to_metadata[sample] = metadata

final_metadata_file.write(first_line)
final_metadata_file.write("\n")
for sample in sample_to_metadata:
    if sample in samples_list:
        #final_metadata_file.write(sample+"\t"+"\t".join(map(str,sample_to_metadata[sample])))
        final_metadata_file.write(sample+"\t"+"\t".join(sample_to_metadata[sample]))
        final_metadata_file.write("\n")

final_metadata_file.close()