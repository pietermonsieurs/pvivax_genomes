#!/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=28 
#SBATCH --time=20:00:00
#SBATCH --job-name=validator

## check the validity of the concatenated vcf file
# cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/filter_noPL/
# /user/antwerpen/205/vsc20587/data/software/vcf-validator/vcf_validator_linux --input PvP01_v1.filtered.no_PL.vcf --report text

cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/filter_core_coveragefilter_noPL/

bcftools view -Ov PvP01_v1.filtered.core.coverage_filter.no_PL.vcf.gz -o PvP01_v1.filtered.core.coverage_filter.no_PL.vcf
/user/antwerpen/205/vsc20587/data/software/vcf-validator/vcf_validator_linux --input PvP01_v1.filtered.core.coverage_filter.no_PL.vcf --report text