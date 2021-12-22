import os,sys

study_id_filename=sys.argv[1]
efetch_metadata_filename=sys.argv[2]
samples_output_filename=sys.argv[3]

study_id_file=open(study_id_filename,'r')
efetch_metadata_file=open(efetch_metadata_filename,'r')
samples_output_file=open(samples_output_filename,'w')

study_list=[]
for line in study_id_file:
    study=line.strip()
    study_list.append(study)

study_id_file.close()
study_id_to_run_id={}
first_line=efetch_metadata_file.readline()
for line in efetch_metadata_file:
    run_accession=line.strip().split("\t")[0]
    study_id=line.strip().split("\t")[20]
    if study_id not in study_id_to_run_id:
        study_id_to_run_id[study_id]=[]
    study_id_to_run_id[study_id].append(run_accession)
print (study_id_to_run_id)
efetch_metadata_file.close()

for study_id in study_id_to_run_id:
    if study_id in study_list:
        #samples_output_file.write("\n".join(study_id_to_run_id[study_id]))
        samples_output_file.write(study_id+"\t"+str(len(study_id_to_run_id[study_id]))+"\t"+",".join(study_id_to_run_id[study_id]))
        samples_output_file.write("\n")