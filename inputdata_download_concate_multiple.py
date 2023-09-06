#!/usr/bin/env python3

import os

my_debug = 0

## file containing all the accession numbers, where the 
## complex accession numbers are separated by a comma
source_file = '/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/fastq_malariagen/MalariaGen_PV_samples.txt'

## specify directories
src_dir = "/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/fastq_malariagen/"
separate_dir = f"{src_dir}/complex/"
concat_dir = f"{src_dir}/concat/"

## create directory
mkdir_command = f"mkdir {separate_dir}"
# os.system(mkdir_command)
mkdir_command2 = f"mkdir {concat_dir}"
# os.system(mkdir_command2)

## open file and skip header
source_fh = open(source_file, 'r', encoding = "ISO-8859-1")
next(source_fh)


for line in source_fh:
    line  = line.rstrip()
    my_debug and print(line)
    line = line.rstrip()
    data = line.split("\t")
    accnr_full = data[8]

    ## some samples do contain multiple accession numbers
    accnr_full = accnr_full.replace('"', '')
    accnr_full = accnr_full.replace(" ", "")
    accnrs = accnr_full.split(",")

    ## check if accession number list is containing more than 
    ## one accession number
    if len(accnrs) > 1:
        print(f"multiple accession numbers found {accnrs}")

        ## create empty concat command 
        cat_R1 = "cat"
        cat_R2 = "cat"

        ## find all the accnrs and move to new location
        for accnr in accnrs: 
            ## create both R1 files
            src_file_R1 = f"{src_dir}/{accnr}_1.fastq.gz"
            src_file_R2 = f"{src_dir}/{accnr}_2.fastq.gz"

            ## move to complex dir
            mv_R1 = f"mv {src_file_R1} {separate_dir}"
            mv_R2 = f"mv {src_file_R2} {separate_dir}"
            # os.system(mv_R1)
            # os.system(mv_R2)

            ## create new files names
            src_file_R1 = f"{separate_dir}/{accnr}_1.fastq.gz"
            src_file_R2 = f"{separate_dir}/{accnr}_2.fastq.gz"            

            ## add to concat command
            cat_R1 = f"{cat_R1} {src_file_R1}"
            cat_R2 = f"{cat_R2} {src_file_R2}"
            
      
        ## set the official name to the first accnr
        out_accnr = accnrs[0]
        out_R1 = f"{concat_dir}/{out_accnr}_1.fastq.gz"
        out_R2 = f"{concat_dir}/{out_accnr}_2.fastq.gz"
        cat_R1 = f"{cat_R1} > {out_R1}"
        cat_R2 = f"{cat_R2} > {out_R2}"
        # print(cat_R1)
        # print(cat_R2)
        os.system(cat_R1)
        os.system(cat_R2)


###### PRINT WARNING ######
print("!!! IMPORTANT !!!")
print("Now move all the fastq_files from the concat directory one directory down. The original fastq files should already been moved with this script to the complex directory. Check whehther the number of R1 fastq files in your source directory corresponds to the number of genomes you plan to add to the database! \n")

