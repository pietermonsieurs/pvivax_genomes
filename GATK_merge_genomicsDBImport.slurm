#!/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=14
#SBATCH --mem-per-cpu=8g
#SBATCH --time=71:00:00
#SBATCH --job-name=GVCFimport


## takes long time, certainly for chromosome 10 until 14. 


## load modules
module load atools
# module load Java
module load GATK

## read in using atools the gvcf files. This can be whatever list of 
## gvcf files you want to read in. Does not have to be the full .gvcf
## list, but can be a subset
# source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_malariagen/gvcf_files.malariagen.csv)
# source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_genomicsdb/chrom_list.csv)
# source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_genomicsdb/chrom_list_crashed.csv)
source <(aenv --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master_merged/chrom_list_regions.csv)

## parameters
output_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master_merged/
threads=14
batch_size=10

chrom=${chrom}

# gatk --java-options "-Xmx50g -Xms50g" \
#     GenomicsDBImport \
#     --genomicsdb-workspace-path ${output_dir}/PvP01_${chrom}_v1.10000.20000/ \
#     --batch-size 50 \
#     -L PvP01_${chrom}_v1:10000-20000 \
#     --sample-name-map mapping_file.PvP01_${chrom}_v1.map \
#     --reader-threads $threads

# gatk --java-options "-Xmx50g -Xms50g" \
#     GenomicsDBImport \
#     --genomicsdb-workspace-path ${output_dir}/PvP01_${chrom}_v1.10000.20000/ \
#     --batch-size 50 \
#     -L PvP01_${chrom}_v1:10000-20000 \
#     --sample-name-map mapping_file.PvP01_${chrom}_v1.test.map \
#     --reader-threads $threads


gatk --java-options "-Xmx100g -Xms100g" \
    GenomicsDBImport \
    --genomicsdb-workspace-path ${output_dir}/PvP01_${chrom}_v1.${start}.${end}/ \
    --batch-size ${batch_size} \
    -L PvP01_${chrom}_v1:${start}-${end} \
    --sample-name-map mapping_file.PvP01_${chrom}_v1.map \
    --reader-threads $threads


## loop over the different chromosomes
# for chrom in {01..14}; do sbatch --export=chrom=${chrom} /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_merge_genomicsDBImport.slurm;done 
# sbatch --export=chrom=MIT /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_merge_genomicsDBImport.slurm
# sbatch --export=chrom=API /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_merge_genomicsDBImport.slurm

## alternative via atools
# cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_genomicsdb
# for chrom in {01..14}; do echo $chrom; done > chrom_list.csv
# vi chrom_list.csv ## add chrom + MIT + API
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_genomicsdb/chrom_list.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_merge_genomicsDBImport.slurm
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_genomicsdb/chrom_list_crashed.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_merge_genomicsDBImport.slurm
# sbatch --array $(arange --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master_merged/chrom_list_regions.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_merge_genomicsDBImport.slurm

