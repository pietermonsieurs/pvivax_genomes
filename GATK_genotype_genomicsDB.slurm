#!/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=4
#SBATCH --time=71:00:00
#SBATCH --job-name=GVCFgenotype


## load modules
module load atools
module load GATK

## read in using atools the gvcf files. This can be whatever list of 
## gvcf files you want to read in. Does not have to be the full .gvcf
## list, but can be a subset
# source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_malariagen/gvcf_files.malariagen.csv)
# source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_genomicsdb/chrom_list.csv)

## read in the list with regions (so not full chromosome, but windows of 100000 nt)
source <(aenv --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master_merged/chrom_list_regions.csv)

## hard coded start and stop
# chrom=01
# start=100001
# end=200000


## parameters
genomicsdb_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master_merged/PvP01_${chrom}_v1.${start}.${end}/
output_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master/
reference_pvivax=/scratch/antwerpen/grp/aitg/arosanas/pvivax_genomes/data/refgenomes/pvivax/PlasmoDB-46_PvivaxP01_Genome.fasta
threads=1


# gatk --java-options "-Xmx16g" GenotypeGVCFs \
gatk --java-options "-Xmx16g" GenotypeGVCFs \
   -R ${reference_pvivax} \
   -V gendb://${genomicsdb_dir} \
   --output ${output_dir}/PvP01_${chrom}_v1.${start}.${end}.vcf.gz \
   --max-alternate-alleles 6 #\
   # -L PvP01_${chrom}_v1:${start}-${end}



## run using atools with chromosomal regions
# module load atools
# sbatch --array $(arange --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master_merged/chrom_list_regions.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_genotype_genomicsDB.slurm


## old implementation: run using atools with full chromosomes
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_genomicsdb/chrom_list.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_genotype_genomicsDB.slurm

 