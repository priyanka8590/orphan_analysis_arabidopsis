import os, sys

expression_filename = sys.argv[1]
low_depth_samples_filename = sys.argv[2]
output_filename = sys.argv[3]

expression_file = open(expression_filename, 'r')
low_depth_samples_file = open(low_depth_samples_filename, 'r')
output_file = open(output_filename, 'w')


sample_to_counts = {}
first_line_exp = expression_file.readline()
for line in expression_file:
    sample = line.strip().split()[0]
    counts = line.strip().split()[1:]
    #print (sample)
    #print (counts)
    sample_to_counts[sample] = counts
    
expression_file.close()
sample_list = []
first_line_samples = low_depth_samples_file.readline()
for line_2 in low_depth_samples_file:
    sample_from_txt = line_2.strip().split()[0]
    print (sample_from_txt)
    #total_size = line_2.strip().split()[1]
    sample_list.append(sample_from_txt)
low_depth_samples_file.close()
#print (sample_list)

output_file.write(first_line_exp)
for sample in sample_to_counts:
    if sample in sample_list:
        output_file.write(sample+"\t"+"\t".join(sample_to_counts[sample]))
        output_file.write("\n")

output_file.close()
        
    