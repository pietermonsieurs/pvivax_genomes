#!/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --job-name=filtervcfs

module load BioTools
module load atools

# source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_clinic/results/ML_hanne/chrom_list.csv)
source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_genomicsdb/chrom_list.csv)


## specify input and output dir, and specify the list containing the 
## samples coming from the Americas, as well as the core genome file
## containing the regions of interest
src_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_combinegvcf/
out_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_combinegvcf/
sample_file=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_combinegvcf/subsets/SAM_ids.csv
core_genome_bed_file=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/Pv_coregenome.bed



## create the vcf_file name. chrom should be read in using
## the atools module
# chrom=API
# chrom=01
chrom=${chrom}

vcf_file_src=${src_dir}/PvP01_${chrom}_v1/PvP01_${chrom}_v1.filtered.vcf.gz
vcf_file_continent=${src_dir}/subsets/PvP01_${chrom}_v1.filtered.SAM.vcf.gz
vcf_file_continent_core=${src_dir}/subsets/PvP01_${chrom}_v1.filtered.SAM.core.vcf.gz


vcf_file_src_snp=${src_dir}/PvP01_${chrom}_v1/PvP01_${chrom}_v1.filtered.snps.vcf.gz
vcf_file_continent_snp=${src_dir}/subsets/PvP01_${chrom}_v1.filtered.SAM.snps.vcf.gz
vcf_file_continent_snp_core=${src_dir}/subsets/PvP01_${chrom}_v1.filtered.SAM.snps.core.vcf.gz




## only select the samples coming from the Americas for all variants
bcftools view -S ${sample_file} \
    ${vcf_file_src} \
    -Oz -o ${vcf_file_continent}

bcftools index -t ${vcf_file_continent}

## only SNPs
bcftools view -S ${sample_file} \
    ${vcf_file_src_snp} \
    -Oz -o ${vcf_file_continent_snp}

bcftools index -t ${vcf_file_continent_snp}

## limit to the core genome for both the files
bcftools view -R ${core_genome_bed_file} \
    ${vcf_file_continent} \
    -Oz -o ${vcf_file_continent_core}
bcftools index ${vcf_file_continent_core}


bcftools view -R ${core_genome_bed_file} \
    ${vcf_file_continent_snp} \
    -Oz -o ${vcf_file_continent_snp_core}
bcftools index ${vcf_file_continent_snp_core}


## create input for module atools 
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_genomicsdb/chrom_list.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_filter_bcftools_percontinent.sh


