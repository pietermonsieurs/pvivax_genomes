#!/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=28
#SBATCH --time=71:00:00
#SBATCH --job-name=vcf_filter

## load modules
module load GATK
module load atools

## read in the chromosomes
source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/chrom_list.csv)
# source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/chrom_list.crashed.csv)


## test chrom
# chrom=API


## specify the input and output directory
vcf_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/
ref_genome=/scratch/antwerpen/grp/aitg/arosanas/pvivax_genomes/data/refgenomes/pvivax/PlasmoDB-46_PvivaxP01_Genome.fasta

## run indexing of the vcf file first
gatk IndexFeatureFile -I ${vcf_dir}/PvP01_${chrom}_v1.vcf.gz

##### SNPs #####
## select SNPs and do basic filtering on the SNPs based
## on the GATK recommended settings

gatk SelectVariants \
-V ${vcf_dir}/PvP01_${chrom}_v1.vcf.gz \
-select-type SNP \
-O ${vcf_dir}/filter/PvP01_${chrom}_v1.snps.vcf.gz

zcat ${vcf_dir}/filter/PvP01_${chrom}_v1.snps.vcf.gz | grep -v "#" | wc


## add filter information
gatk VariantFiltration \
-V ${vcf_dir}/filter/PvP01_${chrom}_v1.snps.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O ${vcf_dir}/filter/PvP01_${chrom}_v1.filter_added.snps.vcf.gz

## select only SNPs that passed the filter - 3310 SNPs selected
gatk SelectVariants \
-R $ref_genome \
-V ${vcf_dir}/filter/PvP01_${chrom}_v1.filter_added.snps.vcf.gz \
--exclude-filtered true \
-O ${vcf_dir}/filter/PvP01_${chrom}_v1.filtered.snps.vcf.gz

zcat ${vcf_dir}/filter/PvP01_${chrom}_v1.filtered.snps.vcf.gz | grep -v "#" | wc


##### INDELS #####
gatk SelectVariants \
-V ${vcf_dir}/PvP01_${chrom}_v1.vcf.gz \
-select-type INDEL \
-O ${vcf_dir}/filter/PvP01_${chrom}_v1.indels.vcf.gz

zcat ${vcf_dir}/filter/PvP01_${chrom}_v1.indels.vcf.gz | grep -v "#" | wc

## add filter information
gatk VariantFiltration \
-V ${vcf_dir}/filter/PvP01_${chrom}_v1.indels.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "FS > 200.0" --filter-name "FS200" \
-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
-O ${vcf_dir}/filter/PvP01_${chrom}_v1.filter_added.indels.vcf.gz

## select only indels that passed the filter - 5347 passed filter
gatk SelectVariants \
-R $ref_genome \
-V ${vcf_dir}/filter/PvP01_${chrom}_v1.filter_added.indels.vcf.gz \
--exclude-filtered true \
-O ${vcf_dir}/filter/PvP01_${chrom}_v1.filtered.indels.vcf.gz

zcat ${vcf_dir}/filter/PvP01_${chrom}_v1.filtered.indels.vcf.gz | grep -v "#" | wc


## merge the SNPs and Indels again to produce one output file
java -Xmx64g -jar ~/data/software/picard/picard.jar   SortVcf I=${vcf_dir}/filter/PvP01_${chrom}_v1.filtered.snps.vcf.gz I=${vcf_dir}/filter/PvP01_${chrom}_v1.filtered.indels.vcf.gz O=${vcf_dir}/filter/PvP01_${chrom}_v1.filtered.vcf.gz
zcat ${vcf_dir}/filter/PvP01_${chrom}_v1.filtered.vcf.gz | grep -v "#" | wc



## run using atools over all chromosomes
# cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/
# echo chrom > chrom_list.csv
# for chrom in {01..14}; do echo $chrom; done >> chrom_list.csv
# echo -e "API\nMIT" >> chrom_list.csv
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/chrom_list.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_filter_vcfs.sh
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/chrom_list.crashed.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_filter_vcfs.sh



