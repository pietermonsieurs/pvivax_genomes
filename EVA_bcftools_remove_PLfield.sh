#!/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=2
#SBATCH --time=20:00:00
#SBATCH --job-name=PLremove

## load modules
module load BioTools
module load atools

## get chromosomes 
source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/chrom_list.csv)


## create input and output files - for the 1532 files, including the 
## ones with too low coverage 
# vcf_in_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/filter/
# vcf_in_file=${vcf_in_dir}/PvP01_${chrom}_v1.filtered.vcf.gz

# vcf_out_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/filter_noPL/
# vcf_out_file=${vcf_out_dir}/PvP01_${chrom}_v1.filtered.no_PL.vcf.gz

## create input and output files - for the 1474 files, with filtering 
## out the bad ones, and corresponding to the mater database
vcf_in_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/filter_core_coveragefilter/
vcf_in_file=${vcf_in_dir}/PvP01_${chrom}_v1.filtered.core.coverage_filter.vcf.gz

vcf_out_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/filter_core_coveragefilter_noPL/
vcf_out_file=${vcf_out_dir}/PvP01_${chrom}_v1.filtered.core.coverage_filter.no_PL.vcf.gz


## convert to a vcf file without PL field
bcftools annotate -x INFO/PL,FORMAT/PL ${vcf_in_file} -Ob -o ${vcf_out_file} --threads 4

## check with the validator script. First you need to unzip the 
## vcf file, as the validator cannot work with gzipped files
vcf_out_file_unzipped=${vcf_out_dir}/PvP01_${chrom}_v1.filtered.no_PL.vcf
bcftools view -Ov -o ${vcf_out_file_unzipped} ${vcf_out_file}
/user/antwerpen/205/vsc20587/data/software/vcf-validator/vcf_validator_linux --input ${vcf_out_file_unzipped} --report text

## run using atools
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/chrom_list.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/EVA_bcftools_remove_PLfield.sh
