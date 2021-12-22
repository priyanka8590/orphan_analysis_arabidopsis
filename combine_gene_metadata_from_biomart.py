import os, sys


def combineBiomartAndOtherMetadata(biomart_file, metadata_file_from_jing, output_file, gene_to_tair_dict, transcript_id_to_cds_start, transcript_id_to_cds_end):
    transcript_id_to_metadata = {}
    transcript_id_list = []
    first_line = biomart_file.readline()
    for line in biomart_file:
        transcript_id, gene_biotype, peptide_seq, gene_name, gene_id, description, chromosome, cds_start, cds_end = line.strip().split("\t")
        transcript_id = transcript_id.strip('"')
        transcript_id_list.append(transcript_id)
        if gene_name == '""':
            gene_name = transcript_id
        #if len(gene_name.split()) == 0:
            #gene_name = transcript_id
        transcript_id_to_metadata[transcript_id] = [gene_biotype, peptide_seq, gene_name, gene_id, description, chromosome, cds_start, cds_end]
        
    first_line_from_file = metadata_file_from_jing.readline()
    qseqid_to_metadata = {}
    for line_2 in metadata_file_from_jing:
        qseqid, strata, chromosome, transcript_start, transcript_end, strand, cds_len, cds_gc, exon_num, mean, median, run, mind, bind, dirinf, maker, braker, araport11 = line_2.strip().split("\t")
        qseqid = qseqid.strip('"')
        qseqid_to_metadata[qseqid] = [strata, chromosome, transcript_start, transcript_end, strand, cds_len, cds_gc, exon_num, mind, bind, dirinf, maker, braker, araport11]
    #print (transcript_id_to_metadata)
    print (qseqid_to_metadata)
    
    """for transcript_id in transcript_id_to_metadata:
        output_file.write(transcript_id+"\t"+"\t".join(transcript_id_to_metadata[transcript_id])+"\t"+"\t".join(qseqid_to_metadata[transcript_id])+"\t"+"\t".join(gene_to_tair_dict[transcript_id]))
        output_file.write("\n")"""
    transcripts_not_in_tair = []
    for transcript_id in transcript_id_to_metadata:
        if transcript_id not in gene_to_tair_dict:
            transcripts_not_in_tair.append(transcript_id)
            
    for transcript_id in transcript_id_to_metadata:
        if transcript_id in gene_to_tair_dict:
            output_file.write(transcript_id+"\t"+"\t".join(transcript_id_to_metadata[transcript_id])+"\t"+"\t".join(qseqid_to_metadata[transcript_id])+"\t"+"\t".join(gene_to_tair_dict[transcript_id])+"\t"+transcript_id_to_cds_start[transcript_id][0]+"\t"+transcript_id_to_cds_end[transcript_id][-1])
            output_file.write("\n")
            
    for transcript in transcripts_not_in_tair:
        output_file.write(transcript+"\t"+"\t".join(transcript_id_to_metadata[transcript])+"\t"+"\t".join(qseqid_to_metadata[transcript])+"\t"+""+"\t"+""+"\t"+""+"\t"+transcript_id_to_cds_start[transcript][0]+"\t"+transcript_id_to_cds_end[transcript][-1])
        output_file.write("\n")
    return qseqid_to_metadata
    
def getMetadataOnlyForEbGenes(qseqid_dict, eb_genes_file, output_file, eb_gene_to_peptide_sequence_dict, transcript_id_to_cds_start, transcript_id_to_cds_end):
    eb_gene_list = []
    first_line = eb_genes_file.readline()
    for line in eb_genes_file:
        eb_gene = line.strip()
        #eb_gene = eb_gene.replace('', '"')
        #eb_gene_quoted = '"'+eb_gene+'"'
        #print (eb_gene_quoted)
        eb_gene_list.append(eb_gene)
    print (qseqid_dict)
    for eb_gene in eb_gene_list:
        output_file.write(eb_gene+ "\t"+"\t".join(qseqid_dict[eb_gene])+"\t"+(eb_gene_to_peptide_sequence_dict[eb_gene]))
        output_file.write("\n")
        
def extractTairDescriptionForGenes(tair_file):
    gene_to_tair_description = {}
    for line in tair_file:
        gene, gene_type, short_description, curator_summary, computational_description = line.strip().split('\t')
        gene = gene.strip('"')
        gene_to_tair_description[gene] = [short_description, curator_summary, computational_description]
    print (gene_to_tair_description)
    return gene_to_tair_description

def createDictOfEbGenesAndPeptideSequence(eb_peptide_sequence_file):
    eb_gene_to_peptide_sequence = {}
    with open(eb_peptide_sequence_file, 'r') as f:
        for line in f:
            if line[0] == ">":
                eb_id = line.strip().split()[0].strip(">")[0:-2]
                #eb_id = '"'+eb_id+'"'
                eb_gene_to_peptide_sequence[eb_id] = f.readline().strip()
            """if line[0:2] == ">AT":
                at_id = line.strip().split()[0].strip(">")[0:-2]
                eb_gene_to_peptide_sequence[at_id] = f.readline().strip()
            else:
                id = line.strip().split()[0].strip(">")[0:]
                eb_gene_to_peptide_sequence[id] = f.readline().strip()"""
    #print (eb_gene_to_peptide_sequence)
    return (eb_gene_to_peptide_sequence)
        
def readGFF3Files(gff3_file):
    transcript_id_to_cds_start = {}
    transcript_id_to_cds_end = {}
    for line_number, line in enumerate(gff3_file):
        if line[0] == "/":continue
        if line[0] == "#":continue
        #print (line, line_number)
        chromosome, source, type, start, end, score, strand, phase, attributes = line.strip().split('\t')
        #print (attributes)
        if type == "CDS":
            transcript_id = attributes.strip().split("=")[1].split(".")[0:2]
            #print (transcript_id)
            transcript_id_final = transcript_id[0]+"."+transcript_id[1]
            #print (transcript_id_final)
            if transcript_id_final not in transcript_id_to_cds_start:
                transcript_id_to_cds_start[transcript_id_final] = []
            transcript_id_to_cds_start[transcript_id_final].append(start)
            if transcript_id_final not in transcript_id_to_cds_end:
                transcript_id_to_cds_end[transcript_id_final] = []
            transcript_id_to_cds_end[transcript_id_final].append(end)
    #print (transcript_id_to_cds_start)
    #print (transcript_id_to_cds_end)
    return transcript_id_to_cds_start, transcript_id_to_cds_end

def mergeGeneMetadataWithCounts(merged_metadata_file, counts_file, output_file):
    gene_to_metadata = {}
    first_line_metadata_file = merged_metadata_file.readline().strip()
    for line in merged_metadata_file:
        gene = line.strip().split("\t")[0]
        metadata = line.strip().split("\t")[1:]
        gene_to_metadata[gene] = metadata
    gene_to_counts = {}
    first_line_counts_file = counts_file.readline().strip()
    for line_2 in counts_file:
        gene = line_2.strip().split("\t")[0]
        counts = line_2.strip().split("\t")[1:]
        gene_to_counts[gene] = counts
    output_file.write(first_line_metadata_file+"\t"+first_line_counts_file)
    output_file.write("\n")
    for gene in gene_to_metadata:
        output_file.write(gene+"\t"+"\t".join(gene_to_metadata[gene])+"\t"+"\t".join(gene_to_counts[gene]))
        output_file.write("\n")
    
    
            

def main():
    biomart_filename = sys.argv[1]
    metadata_filename_from_jing = sys.argv[2]
    output_filename_for_non_eb_genes = sys.argv[3]
    eb_genes_filename = sys.argv[4]
    output_filename_for_eb_genes = sys.argv[5]
    tair_description_filename = sys.argv[6]
    eb_peptide_sequence_filename = sys.argv[7]
    merged_gff3_filename = sys.argv[8]
    merged_gene_metadata_filename = sys.argv[9]
    counts_filename = sys.argv[10]
    merged_gene_metadata_and_counts_output = sys.argv[11]
    
    biomart_file = open(biomart_filename, 'r')
    metadata_file_from_jing = open(metadata_filename_from_jing, 'r')
    output_file = open(output_filename_for_non_eb_genes, 'w')
    eb_genes_file = open(eb_genes_filename, 'r')
    output_file_for_eb_genes = open(output_filename_for_eb_genes, 'w')
    tair_description_file = open(tair_description_filename, 'r')
    gff3_file = open(merged_gff3_filename, 'r')
    merged_gene_metadata_file = open(merged_gene_metadata_filename, 'r')
    counts_file = open(counts_filename, 'r')
    merged_gene_metadata_and_counts_output_file = open(merged_gene_metadata_and_counts_output, 'w')
    transcript_id_to_cds_start, transcript_id_to_cds_end = readGFF3Files(gff3_file)
    #eb_peptide_sequence_file = open(eb_peptide_sequence_filename, 'r')
    eb_gene_to_peptide_dict = createDictOfEbGenesAndPeptideSequence(eb_peptide_sequence_filename)
    gene_to_tair_description = extractTairDescriptionForGenes(tair_description_file)
    qseqid_to_metadata = combineBiomartAndOtherMetadata(biomart_file, metadata_file_from_jing, output_file, gene_to_tair_description, transcript_id_to_cds_start, transcript_id_to_cds_end)
    getMetadataOnlyForEbGenes(qseqid_to_metadata, eb_genes_file, output_file_for_eb_genes, eb_gene_to_peptide_dict, transcript_id_to_cds_start, transcript_id_to_cds_end)
    mergeGeneMetadataWithCounts(merged_gene_metadata_file, counts_file, merged_gene_metadata_and_counts_output_file)
    
if __name__ == "__main__":
    main()