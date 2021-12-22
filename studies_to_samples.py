import os, sys

metadata_filename = sys.argv[1]

metadata_file = open(metadata_filename, 'r')
study_to_run = {}
first_line = metadata_file.readline()
for line in metadata_file:
    study = line.strip().split("\t")[1]
    run = line.strip().split("\t")[0]
    if study not in study_to_run:
        study_to_run[study] = []
    study_to_run[study].append(run)
    
for study in study_to_run:
    if len(study_to_run[study]) == 1:
        #print (study+"\t"+",".join(study_to_run[study]))
        print (study_to_run[study])
    
#print (study_to_run)