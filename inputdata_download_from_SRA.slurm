#!/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --mem-per-cpu=4g
#SBATCH --time=8:00:00
#SBATCH --job-name=SRA_download

module load atools

## read in using atools the gvcf files
# source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/aitg_scratch/arosanas/pvivax_genomes/data/fastq_malariagen/download_url.csv)
source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/fastq_valdivia/sra_numbers_full.fastq_urls.csv)

## move to the correct directory
#cd /user/antwerpen/205/vsc20587/aitg_scratch/arosanas/pvivax_genomes/data/fastq_malariagen/
cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/fastq_valdivia/src_fastq/

wget --tries=0 --wait=120 $url


## run for different fastq-files
## first create a list of URLs by inputting all the accession numbers
## into SRAexplorer. Put "url" as header
# module load atools
# sbatch --array $(arange --data /user/antwerpen/205/vsc20587/aitg_scratch/arosanas/pvivax_genomes/data/fastq_malariagen/download_url.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/inputdata_download_from_SRA.slurm

## alternatively: start from a list of accession numbers and extract
## the download link for the fastq-files using ffq. THis is recommended for
## more complex accession numbers where multiple accession numbers are linked
## to the same sample -> inputdata_download_from_SRA_ffq.sh and inputdata_download_from_SRA_ffq.py
# module load atools
# sbatch --array $(arange --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/fastq_valdivia/sra_numbers_full.fastq_urls.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/inputdata_download_from_SRA.slurm
# tail -n 5 slurm*out | grep "100%" | wc

## alternative
# while read url; do wget $url; done < /user/antwerpen/205/vsc20587/aitg_scratch/arosanas/pvivax_genomes/data/fastq_malariagen/download_url.csv
