#!/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=2
#SBATCH --time=08:00:00
#SBATCH --job-name=filtervcfs

module load BioTools
module load atools

# source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_clinic/results/ML_hanne/chrom_list.csv)
source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/chrom_list.csv)


## specify input and output dir
# src_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_combinegvcf/
# out_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_combinegvcf/

vcf_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/
src_dir=${vcf_dir}/filter_core_coveragefilter_gtmissing/

## select minor allele frequency. In total there are
## ~ 1500 samples --> 0.005 = 0.5% of samples where
## sample should be present = 7.5 samples in database
## of 1500 genomes 
# maf=0.005
# maf=0.010
# maf=0.050
# maf=0.100
# maf=0.200
out_dir=${vcf_dir}/filter_core_coveragefilter_gtmissing_maf${maf}/
mkdir $out_dir



## create the vcf_file name. chrom should be read in using
## the atools module
# chrom=API
# chrom=01
chrom=${chrom}
vcf_file_src=${src_dir}/PvP01_${chrom}_v1.filtered.core.coverage_filter_gtmissing.vcf.gz
vcf_file_filter=${out_dir}/PvP01_${chrom}_v1.filtered.core.coverage_filter_gtmissing.maf${maf}.vcf.gz
vcf_file_filter_snp=${out_dir}/PvP01_${chrom}_v1.filtered.core.coverage_filter_gtmissing.maf${maf}.snps.vcf.gz

# vcf_file_src_snp=${src_dir}/PvP01_${chrom}_v1/PvP01_${chrom}_v1.filtered.snps.vcf.gz
# vcf_file_filter_snp=${src_dir}/PvP01_${chrom}_v1/PvP01_${chrom}_v1.filtered.maf_${maf}.snps.vcf.gz



## SNPs and indels
bcftools view \
    -i '%FILTER="PASS"' \
    -m2 -M2 \
    -q ${maf}:minor \
    ${vcf_file_src} \
    -Oz \
    -o ${vcf_file_filter}

bcftools index -t ${vcf_file_filter}

## filter out the SNPs, and remove indels
bcftools view -v snps -m2 -M2 ${vcf_file_filter} -Oz -o ${vcf_file_filter_snp}
bcftools index -t ${vcf_file_filter_snp}

## only SNPs
# bcftools view \
#     -i '%FILTER="PASS"' \
#     -m2 -M2 \
#     -q ${maf}:minor \
#     ${vcf_file_src_snp} \
#     -Oz \
#     -o ${vcf_file_filter_snp}

# bcftools index -t ${vcf_file_filter_snp}

# check content of the files
for vcf_file in ${out_dir}/*.vcf.gz
do
    # zcat $vcf_file | grep -v "^#" | wc
    zcat $vcf_file | wc
done


## create input for module atools 
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/chrom_list.csv) --export=maf=0.005 /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_filter_bcftools.sh
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/chrom_list.csv) --export=maf=0.010 /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_filter_bcftools.sh
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/chrom_list.csv) --export=maf=0.050 /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_filter_bcftools.sh
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/chrom_list.csv) --export=maf=0.100 /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_filter_bcftools.sh
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/chrom_list.csv) --export=maf=0.200 /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_filter_bcftools.sh