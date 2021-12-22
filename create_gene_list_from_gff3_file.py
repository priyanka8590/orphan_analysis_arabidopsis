import os, sys
import argparse

def parseCommandLineArguments():
    parser = argparse.ArgumentParser(prog="create_gene_lists_from_gff3_file.py",description="")
    
    parser.add_argument("--output_directory","-o",help="Enter the name of the directory where all other operations will be performed",required=True)
    parser.add_argument("--gff3_file","-g",help="Enter the name of the gff3 file",required=True)
    
    return parser.parse_args()

def main():
    commandLineArg=sys.argv
    if len(commandLineArg)==1:
        print("Please use the --help option to get usage information")
    
    options=parseCommandLineArguments()
    createGeneListsFromGff3File(options)
    
def createGeneListsFromGff3File(options):
    gff3_filename = options.gff3_file
    gff3_file = open(gff3_filename, "r")
    species_name = gff3_filename.split("/")[7].split(".")[0]
    print (species_name)
    output_file = options.output_directory+"/"+species_name+"_gene_list.txt"
    fhw = open(output_file, 'w')
    #ID = []
    for line in gff3_file:
        if line[0] == "#":continue
        chromosome, source, type, start, end, score, strand, phase, attributes = line.strip().split("\t")
        if type == "mRNA":
            if "Name" in attributes:
                name = attributes.split(";Name=")[1].split(";")[0].strip()
                fhw.write(name)
                fhw.write("\n")
            elif "ID" in attributes:
                id = attributes.split("ID=")[1].split(";")[0].strip()
                fhw.write(id)
                fhw.write("\n")
            #print (line.strip())
            #ID.append(attributes.split(";")[0].split("=transcript:")[1].strip("'"))
            #ID = attributes.split(";")[0].split("=transcript:")[1].strip("'")
            
            #name = attributes.split(";Name=")[1].split(";")
            #unquoted_ID = ID.strip('"')
            print (name)
            
            
            
            
if __name__ == "__main__":
    main()
            