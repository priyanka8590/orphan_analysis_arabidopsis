import os, sys

filename_from_mind_bind=sys.argv[1]
filename_from_biomart=sys.argv[2]
final_metadata_filename=sys.argv[3]

file_from_mind_bind=open(filename_from_mind_bind, 'r')
file_from_biomart=open(filename_from_biomart,'r')
final_metadata_file=open(final_metadata_filename,'w')

first_line=file_from_mind_bind.readline()
transcipt_to_attributes={}
for line in file_from_mind_bind:
    s_no,qseqid,strata,chromosome,transcript_start,transcript_end,strand,cds_len,cds_gc,exon_num,mean,median,run_2,mind,bind,dirinf,maker,braker,araport11=line.strip().split("\t")
    #qseqid=line.strip().split("\t")[0]
    #print (run_2)
    if qseqid not in transcipt_to_attributes:
        transcipt_to_attributes[qseqid]=[strata,chromosome,transcript_start,transcript_end,strand,cds_len,cds_gc,exon_num,mean,median,run_2,mind,bind,dirinf,maker,braker,araport11]
#print (transcipt_to_attributes["AT1G01100.3"])
    
file_from_mind_bind.close()
transcript_to_attributes_from_biomart={}
for line in file_from_biomart:
    gene_type,peptide,gene_id,description,cds_start,cds_end,transcript_id=line.strip().split("\t")
    #print (cds_end)
    if transcript_id not in transcript_to_attributes_from_biomart:
        transcript_to_attributes_from_biomart[transcript_id]=[gene_type,peptide,gene_id,description,cds_start,cds_end]
#print (transcript_to_attributes_from_biomart)
print (transcript_to_attributes_from_biomart["AT1G01100.3"][3])
file_from_biomart.close()

final_metadata_file.write("Transcript_id"+"\t"+"strata"+"\t"+"chromosome"+"\t"+"transcript_start"+"\t"+"transcript_end"+"\t"+"strand"+"\t"+"cds_len"+"\t"+"cds_gc"+"\t"+"exon_num"+"\t"+"mean"+"\t"+"median"+"\t"+"run_2"+"\t"+"mind"+"\t"+"bind"+"\t"+"dirinf"+"\t"+"maker"+"\t"+"braker"+"\t"+"araport11"+"\t"+"gene_type"+"\t"+"peptide_seq"+"\t"+"gene_id"+"\t"+"description"+"\t"+"cds_start"+"\t"+"cds_end")
final_metadata_file.write("\n")
for transcript_id in transcipt_to_attributes:
    if transcript_id in transcript_to_attributes_from_biomart:
        #print (transcript_id)
        final_metadata_file.write(transcript_id+"\t"+"\t".join(map(str,transcipt_to_attributes[transcript_id]))+"\t"+"\t".join(map(str,transcript_to_attributes_from_biomart[transcript_id])))
        final_metadata_file.write("\n")

final_metadata_file.close()