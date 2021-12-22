import os, sys

transcript_id_filename = sys.argv[1]
output_filename = sys.argv[2]

transcript_id_file = open(transcript_id_filename, 'r')
output_file = open(output_filename, 'w')
first_line = transcript_id_file.readline()

for line in transcript_id_file:
    transcript_id = line.strip()
    gene_id = transcript_id.rsplit(".",1)[0]
    print (transcript_id, gene_id)
    output_file.write(transcript_id+"\t"+gene_id+"\n")
    