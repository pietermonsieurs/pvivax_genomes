#!/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=2
#SBATCH --mem-per-cpu=12g
#SBATCH --time=71:59:00
# #SBATCH --job-name=GATKrun

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=pmonsieurs@itg.be


module load calcua/2020a 
module load atools
module load Java/11
module load GATK/4.1.4.1-Java-8-minimal

# source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/temp/bam_files.csv)
# source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/temp/bam_files.crashed.csv)
source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/temp_malariagen/bam_files.csv)

# threads=2
threads=2

genome_reference=/scratch/antwerpen/grp/aitg/arosanas/pvivax_genomes/data/refgenomes/pvivax/PlasmoDB-46_PvivaxP01_Genome.fasta
# gvcf_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf/
gvcf_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_malariagen/

## test file
# bam_file=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/temp/WGA_2-30.pvivax.bam

bam_file=${bam_file}
sample=${bam_file##*/}
sample=${sample%.pvivax.bam}

gvcf_file=${gvcf_dir}${sample}.g.vcf.gz

# gatk --java-options "-Xmx4g" HaplotypeCaller \
gatk --java-options "-Xmx12g" HaplotypeCaller \
   -R $genome_reference \
   -I $bam_file \
   -O ${gvcf_file} \
   --native-pair-hmm-threads $threads \
   -ERC GVCF



## run using atools
# cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/temp/
# ls $PWD/*.pvivax.bam > bam_files.csv
# vi bam_files.csv # add header bam_file
# echo $PWD/bam_files.csv
# module load atools
# sbatch --array $(arange --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/temp/bam_files.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_run.slurm


## make separate .csv file with the crashed bam-files --> took to
## long due to large size?
# module load atools
# sbatch --array $(arange --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/temp/bam_files.crashed.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_run.slurm


## run on new batch of MalariaGen genomes
# cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/temp_malariagen
# ls $PWD/*.bam  > bam_files.csv
# vi bam_files.csv # add header bam_file
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/temp_malariagen/bam_files.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_run.slurm