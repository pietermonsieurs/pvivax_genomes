#!/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=28 
#SBATCH --time=71:00:00
#SBATCH --job-name=GVCFgenotype

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=pmonsieurs@itg.be


## load modules
module load atools
module load GATK

## read in using atools the gvcf files. This can be whatever list of 
## gvcf files you want to read in. Does not have to be the full .gvcf
## list, but can be a subset
# source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_genomicsdb/chrom_list.csv)

## hard coded chrom instead of via atools
# chrom=06
chrom=11

## parameters
gvcf_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_combinegvcf/PvP01_${chrom}_v1/
reference_pvivax=/scratch/antwerpen/grp/aitg/arosanas/pvivax_genomes/data/refgenomes/pvivax/PlasmoDB-46_PvivaxP01_Genome.fasta
threads=28



if [[ ! -f "${gvcf_dir}/PvP01_${chrom}_v1.g.vcf.gz.tbi" ]]; then
   # echo ${gvcf_dir}/PvP01_${chrom}_v1.g.vcf.gz
   gatk IndexFeatureFile --input ${gvcf_dir}/PvP01_${chrom}_v1.g.vcf.gz
fi


# gatk --java-options "-Xmx16g" GenotypeGVCFs \
gatk --java-options "-Xmx250g -XX:ParallelGCThreads=${threads}" GenotypeGVCFs \
   -R ${reference_pvivax} \
   -V ${gvcf_dir}/PvP01_${chrom}_v1.g.vcf.gz \
   --output ${gvcf_dir}/PvP01_${chrom}_v1.vcf.gz \
   --max-alternate-alleles 2 \
   --max-genotype-count 10 #\
   # -L PvP01_${chrom}_v1:${start}-${end}



## run using atools with full chromosomes
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_genomicsdb/chrom_list.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_genotype_combinegvcfs.slurm

## run using atools with chromosomal regions
# sbatch --array $(arange --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_genomicsdb/chrom_list_regions.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_genotype_genomicsDB.slurm

 