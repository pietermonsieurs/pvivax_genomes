#!/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --mem-per-cpu=4g
#SBATCH --time=20:00:00
#SBATCH --job-name=GVCFsplit

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=pmonsieurs@itg.be


## load modules
module load calcua/2020a 
module load atools
module load Java/11
module load GATK/4.1.4.1-Java-8-minimal

## read in using atools the gvcf files. This can be whatever list of 
## gvcf files you want to read in. Does not have to be the full .gvcf
## list, but can be a subset
# source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master/gvcf_files.csv)
# source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_malariagen/gvcf_files.malariagen.csv)
source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master/gvcf_files.missed.csv)


## parameters
gvcf_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf/
gvcf_perchrom_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_perchrom/
chrom_file=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/refgenomes/pvivax/pvivax_chromosomes.csv

## only do once for creating directories inside the gvcf directory
# while read chrom
# do 
#     chrom_dir=${gvcf_perchrom_dir}/${chrom}/
#     echo $chrom_dir
#     mkdir $chrom_dir
# done < $chrom_file

## test file
# gvcf_file=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf/ERR019040.g.vcf.gz
# gvcf_file=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master/12072018.g.vcf.gz

gvcf_file=${gvcf_file}
while read chrom
do 
    echo "chrom ${chrom}"
    gvcf_file_chrom=${gvcf_file/g.vcf.gz/${chrom}.g.vcf.gz}
    # echo $gvcf_file_chrom
    gvcf_file_chrom=${gvcf_file_chrom/gvcf_master/gvcf_perchrom\/${chrom}}
    echo $gvcf_file_chrom
    gatk SelectVariants -V ${gvcf_file} -L ${chrom} --output ${gvcf_file_chrom}

done < $chrom_file



## run using atools
# cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master/
# ls $PWD/*.g.vcf.gz > gvcf_files.csv
# vi gvcf_files.csv # add header gvcf_file

## copy last lines as maximum number of parameters in atools is 1000
## select last ~ 20 lines and cut 
# tail -28 gvcf_files.csv > gvcf_files.batch2.csv
# vi gvcf_files.batch2.csv ## add header

# echo $PWD/gvcf_files.csv
# echo $PWD/gvcf_files.batch2.csv


# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master/gvcf_files.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_split_gvcf.slurm
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master/gvcf_files.batch2.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_split_gvcf.slurm

## missed files
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master/gvcf_files.missed.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_split_gvcf.slurm


## malariagen gvcf files
# cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_malariagen
# ls $PWD/*.g.vcf.gz > gvcf_files.malariagen.csv
# vi gvcf_files.malariagen.csv ## add header gvcf_file
# echo $PWD/gvcf_files.malariagen.csv
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_malariagen/gvcf_files.malariagen.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_split_gvcf.slurm

## move to the correct directory as output dir was not set correctly
# cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_malariagen/
# out_dir_src=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_perchrom/
# for chrom in {01..14}
# do 
#     # echo $chrom
#     # chrom=API
#     # chrom=MIT
#     out_dir=${out_dir_src}/PvP01_${chrom}_v1/
#     echo $out_dir
#     mv *PvP01_${chrom}* ${out_dir}
# done
