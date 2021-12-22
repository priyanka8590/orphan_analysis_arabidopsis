import os, sys

initial_id_filename = sys.argv[1]
current_id_filename = sys.argv[2]
ids_not_in_initial_filename = sys.argv[3]
initial_id_file = open(initial_id_filename, "r")
current_id_file = open(current_id_filename, "r")
ids_not_in_initial_file = open(ids_not_in_initial_filename, "w")
initial_id_list = []
for line in initial_id_file:
    id = line.strip()
    initial_id_list.append(id)
    
ids_not_in_initial = []
current_id_list = []
for line_2 in current_id_file:
    current_id = line_2.strip()
    current_id_list.append(current_id)
    if current_id not in initial_id_list:
        ids_not_in_initial.append(current_id)
        ids_not_in_initial_file.write(current_id+"\n")
        
