import os,sys

metadata_filename=sys.argv[1]
runs_filename=sys.argv[2]
final_metadata_filename=sys.argv[3]

metadata_file=open(metadata_filename,'r',encoding="ISO-8859-1")
runs_file=open(runs_filename,'r')
final_metadata_file=open(final_metadata_filename,'w')
run_to_metadata={}
first_line=metadata_file.readline()
for line in metadata_file:
    run=line.strip().split("\t")[0]
    #study=line.strip().split("\t")[20]
    metadata=line.strip().split("\t")[1:]
    run_to_metadata[run] = metadata
    
metadata_file.close()
runs_list=[]
for line in runs_file:
    run=line.strip()
    runs_list.append(run)
    
runs_file.close()
final_metadata_file.write(first_line)
final_metadata_file.write("\n")
for run in run_to_metadata:
    if run in runs_list:
        final_metadata_file.write(run+"\t"+"\t".join(map(str,run_to_metadata[run])))
        final_metadata_file.write("\n")
final_metadata_file.close()
    
