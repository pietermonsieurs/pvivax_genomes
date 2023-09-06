#!/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=28
#SBATCH --time=70:00:00
#SBATCH --job-name=bwa_run



## load necessary module ##
module load BioTools
module load BWA
module load Java

## set parameters
picard_jar=/user/antwerpen/205/vsc20587/data/software/picard/picard.jar
picard_bin="java -jar ${picard_jar}"
reference_human=/scratch/antwerpen/grp/aitg/arosanas/pvivax_genomes/data/refgenomes/human/GRCh38_latest_genomic.fasta
reference_pvivax=/scratch/antwerpen/grp/aitg/arosanas/pvivax_genomes/data/refgenomes/pvivax/PlasmoDB-46_PvivaxP01_Genome.fasta
# bwa_dir=/scratch/antwerpen/grp/aitg/arosanas/pvivax_genomes/results/temp/
# bwa_dir=/scratch/antwerpen/grp/aitg/arosanas/pvivax_genomes/results/temp_malariagen/
bwa_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/temp_valdivia/
threads=28


## fastq_file_R1: should be read in via the parameter setting but below
## an example fastq_file_R1 file
# fastq_file_1=/scratch/antwerpen/grp/aitg/arosanas/pvivax_genomes/data/fastq_new/ERR5740834_1.fastq.gz
fastq_file_1=${fastq_file_1}
echo "fastq_file_1..."
echo $fastq_file_1

## create the fastq_file_2 and bam file
file_prefix_full=${fastq_file_1%_1.fastq.gz}
fastq_file_2=${file_prefix_full}_2.fastq.gz
sample=${file_prefix_full##*/}
bam_file_prefix=${bwa_dir}/${sample}

echo $fastq_file_1
echo $fastq_file_2
echo $sample
echo $bam_file_prefix


## align the reads first to the human genome to filter out the
## native human reads
bwa mem \
	-R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" \
	-t $threads \
	$reference_human \
	$fastq_file_1 \
	$fastq_file_2 \
	| samtools view -@ $threads -b > ${bam_file_prefix}.human.bam

## run flagstat for statistics
samtools flagstat -@ ${threads} ${bam_file_prefix}.human.bam > ${bam_file_prefix}.human.flagstat

## filter out the proper paired reads - only keep the non-mapped reads. Afterwards
## sort again this bam file with unmapped reads
samtools view -@ ${threads} -b -F 2 -o ${bam_file_prefix}.no_human.bam ${bam_file_prefix}.human.bam
samtools sort -@ ${threads} -n -o ${bam_file_prefix}.no_human.sorted.bam ${bam_file_prefix}.no_human.bam

## extract the fastq files from the unmapped bam file
samtools fastq -@ ${threads} -N ${bam_file_prefix}.no_human.sorted.bam -1 ${bam_file_prefix}.no_human_1.fastq.gz -2 ${bam_file_prefix}.no_human_2.fastq.gz


## map the unmapped reads versus the Pvivax genome
bwa mem \
	-R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" \
	-t $threads \
	$reference_pvivax \
	${bam_file_prefix}.no_human_1.fastq.gz \
	${bam_file_prefix}.no_human_2.fastq.gz \
	| samtools sort -@ $threads -o ${bam_file_prefix}.pvivax.bam

## flagastat and index on pvivax bam
samtools flagstat -@ $threads ${bam_file_prefix}.pvivax.bam > ${bam_file_prefix}.pvivax.flagstat
samtools index -@ $threads ${bam_file_prefix}.pvivax.bam


## remove duplicates from the .pvivax.bam file
## set parameters
$picard_bin MarkDuplicates REMOVE_DUPLICATES=true I=${bam_file_prefix}.pvivax.bam O=${bam_file_prefix}.pvivax.removedups.bam M=${bam_file_prefix}.pvivax.markdups.txt
samtools index ${bam_file_prefix}.pvivax.removedups.bam
samtools flagstat ${bam_file_prefix}.pvivax.removedups.bam > ${bam_file_prefix}.pvivax.removedups.flagstat


## clean up 
rm ${bam_file_prefix}.human.bam
rm ${bam_file_prefix}.no_human.bam
rm ${bam_file_prefix}.pvivax.bam
rm ${bam_file_prefix}.no_human.sorted.bam
rm ${bam_file_prefix}.no_human_1.fastq.gz
rm ${bam_file_prefix}.no_human_2.fastq.gz


## run script per fastq file\
# cd /scratch/antwerpen/grp/aitg/arosanas/pvivax_genomes/data/fastq_new
# sbatch --export=fastq_file_1=/scratch/antwerpen/grp/aitg/arosanas/pvivax_genomes/data/fastq_new/ERR5740834_1.fastq.gz /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/BWA_run.sh

# for fastq_file_1 in ${PWD}/*_1.fastq.gz; do sbatch --export=fastq_file_1=${fastq_file_1} /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/BWA_run.sh; done

## run through the crashed file
# while read fastq_file_1; do echo ${fastq_file_1}; done < /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/fastq_new/crashed_fastqfiles.csv
# while read fastq_file_1; do sbatch --export=fastq_file_1=${fastq_file_1} /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/BWA_run.sh; done < /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/fastq_new/crashed_fastqfiles.csv
# while read fastq_file_1; do sbatch --export=fastq_file_1=${fastq_file_1} /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/BWA_run.sh; done < /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/fastq_new/crashed_fastqfiles.run2.csv


## run for genomes of Erin
# cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/fastq_erin
# for fastq_file_1 in ${PWD}/*_1.fastq.gz; do sbatch --export=fastq_file_1=${fastq_file_1} /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/BWA_run.sh; done


## run for genomes of Marcelo
# cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/fastq_marcelo
# for fastq_file_1 in ${PWD}/*_1.fastq.gz; do sbatch --export=fastq_file_1=${fastq_file_1} /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/BWA_run.sh; done

## run for genomes of MalariaGen data 
# cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/fastq_malariagen
# for fastq_file_1 in ${PWD}/*_1.fastq.gz; do sbatch --export=fastq_file_1=${fastq_file_1} /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/BWA_run.sh; done


## run for genomes of Hugo Valdivia paper (data Peru 2022) 
# cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/fastq_valdivia/
# for fastq_file_1 in ${PWD}/*_1.fastq.gz; do sbatch --export=fastq_file_1=${fastq_file_1} /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/BWA_run.sh; done