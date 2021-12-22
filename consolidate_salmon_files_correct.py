import os, glob, sys

run_filename = sys.argv[1]
raw_count_converged_filename = sys.argv[2]
tpm_converged_filename = sys.argv[3]
run_file = open(run_filename, "r")
raw_count_converged_file = open(raw_count_converged_filename, "w")
tpm_converged_file = open(tpm_converged_filename, "w")
gene_to_raw_count = {}
gene_to_tpm = {}
run_list = []
for line in run_file:
    run = line.strip()
    run_list.append(run)
    salmon_quant_file = open("/work/LAS/mash-lab/bhandary/bind_prediction_expression/araport11_salmon_quant_4786/"+run+"_salmon/quant.genes.sf")
    first_line = salmon_quant_file.readline()
    for line_1 in salmon_quant_file:
        #if "TPM" in line:continue
        gene = line_1.strip().split("\t")[0]
        tpm = line_1.strip().split("\t")[3]
        raw_count = line_1.strip().split("\t")[4]
        print (tpm, raw_count)
        if gene not in gene_to_tpm:
            gene_to_tpm[gene] = []
        gene_to_tpm[gene].append(tpm)
        if gene not in gene_to_raw_count:
            gene_to_raw_count[gene] = []
        gene_to_raw_count[gene].append(raw_count)


raw_count_converged_file.write("Gene"+"\t"+"\t".join(run_list)+"\n")
for gene in gene_to_raw_count:
    raw_count_converged_file.write(gene+"\t"+"\t".join(gene_to_raw_count[gene])+"\n")

tpm_converged_file.write("Gene"+"\t"+"\t".join(run_list)+"\n")
for gene in gene_to_tpm:
    tpm_converged_file.write(gene+"\t"+"\t".join(gene_to_tpm[gene])+"\n")

run_file.close()
raw_count_converged_file.close()
tpm_converged_file.close()