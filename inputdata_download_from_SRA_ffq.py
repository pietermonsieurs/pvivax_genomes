#!/usr/bin/env python3

import json

my_debug = 0
accession_number_file = '/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/fastq_valdivia/sra_numbers_full.csv'
accession_number_file_out = '/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/fastq_valdivia/sra_numbers_full.ffq_SRA.csv'
fastq_file_url = '/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/fastq_valdivia/sra_numbers_full.fastq_urls.csv'
data_dir = '/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/fastq_valdivia/'

## open input file with accession numbers, and adapt
## header for new output
ac_fh = open(accession_number_file, 'r')
header = ac_fh.readline()
header = header.rstrip()
header_new = f"{header}\texp_number\trun_number\tfastq_file_1\tmd5_file_1\tfastq_file_2\tmd5_file_2\n"

## open output file which will have the same information
## as the input files, but some lines duplicated depending
## on number of experiments + columns added to specify the 
## URL and md5 values
out_fh = open(accession_number_file_out, 'w')
out_fh.write(header_new)

## open files to write fastq-files to 
fastq_fh = open(fastq_file_url, 'w')
header = "url\n"
fastq_fh.write(header)


multiple_runs  = 0

for line in ac_fh:
    
    ## extract the accession number from the input file
    ## and create the .json file
    line = line.rstrip()
    data = line.split("\t")
    accnr = data[2]
    json_file = f"{data_dir}/{accnr}.json"
    print(f"json file {json_file}")

    ## load the json file and scroll through it
    json_fh = open(json_file, 'r')
    json_data = json.load(json_fh)
    # print(len(json_data[accnr]['experiments']))
    
    
    if len(json_data[accnr]['experiments']) == 1:
        print(list(json_data[accnr]['experiments'].keys())[0])
        
        ## extract the accession number of the experiment
        sra_number = list(json_data[accnr]['experiments'].keys())[0]

        ## add sra_number to additional columns
        columns_to_add = f"{sra_number}"

        ## extract the accession number of the runs. It assumes
        ## that every experiment only has one run in the json file, whic
        ## is the case for the samples of Valdivia et al. 
        # if len(list(json_data[accnr]['experiments'][sra_number]["runs"].keys())) > 1:
        #     print("MULTIPLE RUNS!!")
        run_number = list(json_data[accnr]['experiments'][sra_number]["runs"].keys())[0]
        columns_to_add = f"{columns_to_add}\t{run_number}"

        for fastq_file_info in json_data[accnr]['experiments'][sra_number]["runs"][run_number]["files"]["ftp"]:
            my_debug and print(fastq_file_info['url'])
            my_debug and print(fastq_file_info['md5'])
            columns_to_add = f"{columns_to_add}\t{fastq_file_info['url']}\t{fastq_file_info['md5']}"
            
            ## write also the fastq file out
            fastq_line_out = f"{fastq_file_info['url']}\n"
            fastq_fh.write(fastq_line_out)
        
        print(columns_to_add)
        line_out = f"{line}\t{columns_to_add}\n"
        out_fh.write(line_out)

            

    elif len(json_data[accnr]['experiments']) > 1:
        multiple_runs = multiple_runs + 1
        for exp_list in json_data[accnr]['experiments']:

            # print(exp_list)
            sra_number = list(exp_list.keys())[0]
            print(sra_number)
            
            ## add sra_number to additional columns
            columns_to_add = f"{sra_number}"

            ## extract the accession number of the runs. It assumes
            ## that every experiment only has one run in the json file, whic
            ## is the case for the samples of Valdivia et al. 
            run_number = list(exp_list[sra_number]['runs'].keys())[0]
            columns_to_add = f"{columns_to_add}\t{run_number}"

            for fastq_file_info in exp_list[sra_number]['runs'][run_number]["files"]["ftp"]:
                my_debug and print(fastq_file_info['url'])
                my_debug and print(fastq_file_info['md5'])
                columns_to_add = f"{columns_to_add}\t{fastq_file_info['url']}\t{fastq_file_info['md5']}"

                ## write also the fastq file out
                fastq_line_out = f"{fastq_file_info['url']}\n"
                fastq_fh.write(fastq_line_out)


            print(columns_to_add)
            line_out = f"{line}\t{columns_to_add}\n"
            out_fh.write(line_out)

    # print(json_data[accnr])


print(f"multiple runs -> {multiple_runs}")
