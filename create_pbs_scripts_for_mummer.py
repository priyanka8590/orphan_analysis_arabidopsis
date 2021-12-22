import os, sys

target_species = ["Arabidopsis_lyrata", "Brassica_napus", "Brassica_olerecea", "An-1", "Cvi", "C24", "Eri1", "Kyo", "Ler", "Sha", "Eutrema_salsugineum"]
new_file = "/work/LAS/mash-lab/bhandary/analysis_orphan_gene_aim_3/fagin_analysis/pbs_scripts/script_for_running_mummer.pbs"
fhw=open(new_file,"w")

fhw.write("""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=75:00:00
#SBATCH --mem=90G
#SBATCH --constraint=AVX2
#SBATCH --partition=biocrunch
#SBATCH --mail-user=bhandary@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --error=/work/LAS/mash-lab/bhandary/analysis_orphan_gene_aim_3/script_for_running_mummer.error
#SBATCH --output=/work/LAS/mash-lab/bhandary/analysis_orphan_gene_aim_3/script_for_running_mummer.output""")
fhw.write("\n")

for i in target_species:
    fhw.write("query=Arabidopsis_thaliana")
    fhw.write("\n")
    fhw.write("target="+str(i))
    fhw.write("\n")
    fhw.write("prefix=${query}_${target}")
    fhw.write("\n")
    fhw.write("nucmer --mum -t 36 -p ${prefix} /work/LAS/mash-lab/bhandary/analysis_orphan_gene_aim_3/fagin_analysis/${query}.fna /work/LAS/mash-lab/bhandary/analysis_orphan_gene_aim_3/fagin_analysis/${target}.fna")
    fhw.write("\n")
    fhw.write("show-coords -T -r -c -l ${prefix}.delta > ${prefix}.coords")
    fhw.write("\n")