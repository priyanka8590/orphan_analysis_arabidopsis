import os, sys

samples_filename = sys.argv[1]
transposed_count_filename = sys.argv[2]
subsetted_filename = sys.argv[3]
samples_file = open(samples_filename, "r")
count_file = open(transposed_count_filename, "r")
subsetted_file = open(subsetted_filename, "w")

sample_list = []
for line in samples_file:
    sample = line.strip()
    sample_list.append(sample)
    
first_line = count_file.readline()
subsetted_file.write(first_line)
#subsetted_file.write("\n")
sample_to_count = {}
for line_2 in count_file:
    sample = line_2.strip().split()[0]
    counts = line_2.strip().split()[1:]
    sample_to_count[sample] = counts
    """if sample not in sample_to_count:
        sample_to_count[sample]=[]
    sample_to_count[sample].append(counts)"""

#print (sample_to_count)
#print (sample_to_count["SRR11637987"])
for sample in sample_to_count:
    if sample in sample_list:
        subsetted_file.write(sample+"\t"+"\t".join(map(str,sample_to_count[sample])))
        subsetted_file.write("\n")