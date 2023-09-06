#!/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=28 
#SBATCH --time=71:00:00
#SBATCH --job-name=GVCFcomb

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=pmonsieurs@itg.be


## load modules
module load atools
module load GATK

## read in using atools the gvcf files. This can be whatever list of 
## gvcf files you want to read in. Does not have to be the full .gvcf
## list, but can be a subset
# source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_malariagen/gvcf_files.malariagen.csv)
# source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_genomicsdb/chrom_list.csv)

## read in the list with regions (so not full chromosome, but windows of 100000 nt)
# source <(aenv --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_genomicsdb/chrom_list_regions.csv)

## run for one specific chromosome instead of via atools
# chrom=06
chrom=11

## parameters
gvcf_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_perchrom/PvP01_${chrom}_v1/
out_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_combinegvcf/PvP01_${chrom}_v1/
reference_pvivax=/scratch/antwerpen/grp/aitg/arosanas/pvivax_genomes/data/refgenomes/pvivax/PlasmoDB-46_PvivaxP01_Genome.fasta
threads=28

## create output directory
mkdir ${out_dir}

## create a list of all gvcf files
cd $gvcf_dir
ls $PWD/*.g.vcf.gz > vcf_files.list


# gatk --java-options "-Xmx16g" GenotypeGVCFs \
gatk --java-options "-Xmx250G -XX:ParallelGCThreads=${threads}" CombineGVCFs \
   -R ${reference_pvivax} \
   --variant ${gvcf_dir}/vcf_files.list \
   --output ${out_dir}/PvP01_${chrom}_v1.g.vcf.gz \
   --exclude-intervals PvP01_${chrom}_v1:1124300-1124800 
   # --exclude-intervals PvP01_${chrom}_v1:107470-107800


## run using atools with full chromosomes
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_genomicsdb/chrom_list.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_combinegvcfs.slurm

## run using atools with chromosomal regions
# sbatch --array $(arange --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_genomicsdb/chrom_list_regions.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_genotype_genomicsDB.slurm

 