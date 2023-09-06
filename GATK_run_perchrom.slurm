#!/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=2
#SBATCH --mem-per-cpu=12g
#SBATCH --time=71:59:00
#SBATCH --job-name=GATKrun

module load atools
module load Java
module load GATK

# source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/bam_master/bam_files.csv)
# source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/bam_master/bam_files_test.csv)
# source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/bam_master/bam_files.batch1.csv)
source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/bam_master/bam_files_missed.csv)

echo "BAM FILE 1 ${bam_file}"

threads=10

genome_reference=/scratch/antwerpen/grp/aitg/arosanas/pvivax_genomes/data/refgenomes/pvivax/PlasmoDB-46_PvivaxP01_Genome.fasta
gvcf_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master/


## test file
# bam_file=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/bam_master/ERR054088.bam

bam_file=${bam_file}
sample=${bam_file##*/}
# sample=${sample%.pvivax.bam}
sample=${sample%.bam}

echo "BAM FILE 2 ${bam_file}"


## first create a chromosome list by looping over the different 
## chromosomes
chrom_list=(01)
for chrom in {02..14}
do
    # echo $chrom
    chrom_list+=(${chrom})
done

chrom_list+=("API")
chrom_list+=("MIT")
echo ${chrom_list[@]}

## if first run, first create directory per chromosome
## in the gvcf_master directory
# cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master
# for chrom in ${chrom_list[@]}; do mkdir PvP01_${chrom}_v1; done

## run over the chromosome list, and run the
## HaplotypeCaller per chrom
echo "BAM FILE 3 ${bam_file}"

for chrom in ${chrom_list[@]}
do
   chrom_full=PvP01_${chrom}_v1
   gvcf_file=${gvcf_dir}/${chrom_full}/${sample}.${chrom_full}.g.vcf.gz

   # gatk --java-options "-Xmx4g" HaplotypeCaller \
   gatk --java-options "-Xmx4g" HaplotypeCaller \
      -R $genome_reference \
      -I $bam_file \
      --intervals $chrom_full \
      -O ${gvcf_file} \
      --native-pair-hmm-threads $threads \
      -ERC GVCF
done



## run using atools
# cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/bam_master/
# ls $PWD/*.bam > bam_files.csv
# vi bam_files.csv # add header bam_file
# echo $PWD/bam_files.csv
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/bam_master/bam_files.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_run_perchrom.slurm
# sbatch --array $(arange --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/bam_master/bam_files_test.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_run_perchrom.slurm
# sbatch --array $(arange --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/bam_master/bam_files.batch1.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_run_perchrom.slurm

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


## run on samples with missing gvcf file after checking the 
## integrity via the GATK_run_perchrom_checkintegrity.sh
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/bam_master/bam_files_missed.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/GATK_run_perchrom.slurm   