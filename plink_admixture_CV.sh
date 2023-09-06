#!/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=7
#SBATCH --time=71:00:00


## read in the K-value using the atools module
## test value: K=2
module load atools
source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/plink/K_values.csv)


## settings for the software to be used
admixture=/user/antwerpen/205/vsc20587/data/software/admixture_linux-1.3.0/admixture
cross_validation=10

## input parameters of the .vcf file and output 
threads=7
maf=0.005
# out_prefix=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/plink/PvP01_v1.filtered.core.coverage_filter_gtmissing.maf0.005.snps.ldpruned.filtered 
out_prefix=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/plink/PvP01_v1.filtered.core.coverage_filter_gtmissing.maf0.005.snps.SAM.ldpruned.filtered


$admixture -j${threads} --cv=${cross_validation} ${out_prefix}.bed $K  > ${out_prefix}.K${K}.log


## create file with differnet values of K
# cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/plink
# for K in {2..50}; do echo $K; done > K_values.csv
# for K in {22..50}; do echo $K; done > K_values.csv
# echo $PWD/K_values.csv
# vi K_values.csv
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/plink/K_values.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/plink_admixture_CV.sh

## afterwards concatenate the different CV error values



