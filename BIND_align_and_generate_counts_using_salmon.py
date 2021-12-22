import os

def main():
    all_samples=open("/work/LAS/mash-lab/bhandary/bind_prediction_expression/4786_srr_ids","r").read().split("\n")
    i=0

    cmd="salmon "
    cmd+=" index -t /work/LAS/mash-lab/bhandary/bind_prediction_expression/salmon_quant/BIND.gentrome.fa "
    cmd+=" -d /work/LAS/mash-lab/bhandary/bind_prediction_expression/salmon_quant/decoys.txt "
    cmd+=" -p 64 "
    cmd+=" -i /work/LAS/mash-lab/bhandary/bind_prediction_expression/salmon_quant/BIND_salmon_index "
    cmd+=" --keepDuplicates "
    #print(cmd)
    os.system(cmd)

    #return
    while True:
        if i>len(all_samples):break
        open("/work/LAS/mash-lab/bhandary/bind_prediction_expression/4786_srr_ids_temp","w").write("\n".join(all_samples[i:i+5]))

        for Run in all_samples[i:i+5]:

            cmd="salmon quant "
            cmd+=" -i /work/LAS/mash-lab/bhandary/bind_prediction_expression/salmon_quant/BIND_salmon_index "
            cmd+=" -l A "
            cmd+=" -p 64 "
            cmd+=" --validateMappings "
            cmd+=" --gcBias "
            cmd+=" -g /work/LAS/mash-lab/bhandary/analysis_orphan_gene_aim_3/transcript_to_gene_mapping "
            cmd+=" --noBiasLengthThreshold "
            cmd+=" -1 /work/LAS/mash-lab/bhandary/bind_prediction_expression/newest_samples_folder/"+Run+"_1.fastq "
            cmd+=" -2 /work/LAS/mash-lab/bhandary/bind_prediction_expression/newest_samples_folder/"+Run+"_2.fastq "
            cmd+=" -o /work/LAS/mash-lab/bhandary/bind_prediction_expression/BIND_salmon_quant_4786/"+Run+"_salmon "
            cmd+=" 2> "+"/work/LAS/mash-lab/bhandary/bind_prediction_expression/BIND_salmon_quant_4786/"+Run+"_salmon.error "
            #print (cmd)
            os.system(cmd)

        #os.system("rm -rf /work/LAS/mash-lab/bhandary/old_analysis_regulon_prediction/open_reading_frame/"+Run+"_STAR*")
        #os.system("rm "+"/work/LAS/mash-lab/bhandary/old_analysis_regulon_prediction/open_reading_frame/"+Run+"_*fastq ")
        i+=5

if __name__ == "__main__":
    main()