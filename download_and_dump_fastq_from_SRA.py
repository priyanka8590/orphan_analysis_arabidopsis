#!/usr/bin/env python

"""
Parallely downloads and converts NCBI .sra files to FASTQ
"""

import argparse
import os
import sys
#from pathlib import Path
import multiprocessing

def parseCommandLineArguments():
    """
    Parses the arguments provided through command line.
    Launch python download_and_dump_fastq_from_SRA.py --help for more details
    """
    parser = argparse.ArgumentParser(prog="download_and_dump_fastq_from_SRA.py",description="Parallel download of fastq data from NCBI.")
    parser.add_argument("--sra","-s",help="Please enter the name of the file which has all the SRAs listed one per line",required=True)
    parser.add_argument("--output","-o",help="Please enter the name of the output directory.",required=True)
    parser.add_argument("--cpu","-n",help="Enter the number of CPUs to be used.",default=1)
    return parser.parse_args()

def readSRAfilesToBeDownloaded(filename):
    """
    Reads and returns a list of the SRA ids to be downloaded
    """
    return list(set([name.strip() for name in open(filename,"r").read().split("\n")]))

def downloadSRAFile(allinput):
    sra,default_path_to_download,output_directory=allinput
    os.system("prefetch -X 104857600 -O "+output_directory+"/"+" "+sra+" 2> "+output_directory+"/"+sra+".error")
    #print("prefetch -O "+output_directory+"/"+" "+sra+" 2> "+output_directory+"/"+sra+".error")
    #os.system("mv "+default_path_to_download+"/"+sra+".sra "+output_directory+"/")
    cmd="fastq-dump -X 1 -Z  --split-spot "+output_directory+"/"+sra+".sra|wc -l > "+output_directory+"/"+sra+".temp"
    os.system(cmd)
    if int(open(output_directory+"/"+sra+".temp").read())==4:
        pair="single"
    else:
        pair="paired"
    cmd="fastq-dump --defline-seq '@$sn[_$rn]/$ri' --outdir "+output_directory+" --split-files "+output_directory+"/"+sra+".sra"
    print (pair)
    os.system(cmd)
    if pair=="single":
        os.system("mv "+output_directory+"/"+sra+"_1.fastq "+output_directory+"/"+sra+".fastq ")
    os.system("rm "+output_directory+"/"+sra+".sra "+output_directory+"/"+sra+".temp")
    
def downloadSRAFilesAndConvertToFastq(SRAs,default_path_to_download,n,output_directory):
    """
    Downloads the sra files and converts to fastq
    """
    cmd="mkdir -p "+output_directory
    os.system(cmd)
    pool = multiprocessing.Pool(processes=int(n))
    allinputs=[]
    os.system("rm -rf "+output_directory+"/*lock")
    os.system("rm -rf "+output_directory+"/*tmp")
    os.system("rm -rf "+output_directory+"/*error")
    os.system("rm -rf "+output_directory+"/*temp")
    for sra in SRAs:
        if os.path.exists(output_directory+"/"+sra+".fastq")==True or (os.path.exists(output_directory+"/"+sra+"_1.fastq")==True and os.path.exists(output_directory+"/"+sra+"_2.fastq")==True):
            if os.path.exists(output_directory+"/"+sra+"_1.fastq")==True and os.path.exists(output_directory+"/"+sra+"_2.fastq")==False:
                os.system("mv "+output_directory+"/"+sra+"_1.fastq "+output_directory+"/"+sra+".fastq")
            continue
        allinputs.append([sra,default_path_to_download,output_directory])
    pool.map(downloadSRAFile,allinputs)
    
def verifyOutput(output_directory,SRAs):
    """
    Verify the downloads
    """
    for sra in SRAs:
        if os.path.exists(output_directory+"/"+sra+".fastq")==True:
            continue
        elif os.path.exists(output_directory+"/"+sra+"_1.fastq")==True and os.path.exists(output_directory+"/"+sra+"_2.fastq")==True:
            continue
        else:
            return 1
        return 0
          
def main():
    commandLineArg=sys.argv
    if len(commandLineArg)==1:
        print("Please use the --help option to get usage information")
    options=parseCommandLineArguments()
    SRAs=readSRAfilesToBeDownloaded(options.sra)
    new_SRAs=[s for s in SRAs if s!=""]
    SRAs=new_SRAs
    #print(SRAs)
    #home = str(Path.home())
    #home="/home/sagnik.banerjee"
    #home="/N/dc2/scratch/prbhan"
    #default_path_to_download=home+"/ncbi/public/sra/"
    default_path_to_download=""
    while verifyOutput(options.output,SRAs)==1:
        downloadSRAFilesAndConvertToFastq(SRAs,default_path_to_download,int(options.cpu),options.output)
    

if __name__ == "__main__":
    main()
