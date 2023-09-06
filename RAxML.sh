#!/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=28 
#SBATCH --time=20:00:00
#SBATCH --job-name=RAxML

## set the number of threads depending on the amount of lrocs
## specified above
threads=16

## load both modules needed for MUSCLE and RAxML
# module load MUSCLE
module load RAxML-NG 

## define the fasta file, either by going to the directory
## and selecting a file, or by the worker module (comment those lines)
# cd /user/antwerpen/205/vsc20587/scratch/ext_cmpg_primer/results/raxml
# fasta_file=A2gamF_A3_gamR_Tambadou2014.amplicons.fasta
# fasta_file=A3_A7_AyusoSacido2005.amplicons.fasta
# fasta_file=ApIF_TpIR_Tapi2010.amplicons.fasta
# fasta_file=MTR_MTF_Tambadou2014.amplicons.fasta

## run MUSCLE
# muscle -in $fasta_file -phyiout ${fasta_file}.phy

## path for the conversion tool of vcf to Phylip
vcf2phylip=/user/antwerpen/205/vsc20587/data/software/vcf2phylip/vcf2phylip.py

## input data and convert to phylip: WORLD
# vcf_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_combinegvcf/test_plink/
# vcf_prefix=PvP01_v1.filtered.maf_0.200.snps.core.plink_filtered

## input data and convert to phylip: SAM
vcf_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/filter_core_coveragefilter_gtmissing_maf0.005
vcf_prefix=PvP01_v1.filtered.core.coverage_filter_gtmissing.maf0.005.snps.SAM


file_prefix=${vcf_dir}/${vcf_prefix}

$vcf2phylip  -i ${file_prefix}.vcf.gz --output-folder $vcf_dir --output-prefix ${file_prefix}.phy

## output file has sometimes different naming, then make symbolic
## linkt to assumed file name
if [ ! -f "${file_prefix}.phy" ]; then
    echo "${file_prefix}.phy does not exist."
    ln -s ${file_prefix}.phy.min4.phy ${file_prefix}.phy
fi

## RAxML - check the input data
# raxml-ng --check --msa ${fasta_file}.phy --model LG+G8+F --prefix ${fasta_file}
raxml-ng --check --msa ${file_prefix}.phy --model GTR+G --prefix ${file_prefix}.T1

if [ ! -f "${file_prefix}.T1.raxml.reduced.phy" ]; then
    echo "${file_prefix}.T1.raxml.reduced.phy does not exist."
    ln -s ${file_prefix}.phy ${file_prefix}.T1.raxml.reduced.phy
fi

## RAxML - convert to binary format and do estimates of CPU
raxml-ng --parse --msa ${file_prefix}.T1.raxml.reduced.phy --model GTR+G --prefix ${file_prefix}.T2


## RAxML run the program
# raxml-ng --msa ${fasta_file}.T2.raxml.rba --model GTR+G --prefix ${fasta_file}.T3 --threads $threads --seed 2 
raxml-ng --search1 --msa ${file_prefix}.T2.raxml.rba --model GTR+G --prefix ${file_prefix}.T3 --threads $threads --seed 2 


## RAxML, run with bootstrapping
# raxml -f a -# 20 -m PROTGAMMAAUTO -p 12345 -x 12345 -s toy_dataset_PSII_protein_aligned.phy -n toy_dataset_PSII_protein.tree
# raxml-ng --all --msa ${fasta_file}.phy --model LG+G8+F --tree pars{10} --bs-trees 200 --threads $threads
# raxml-ng --all --msa ${fasta_file}.phy --model LG+G8+F --tree pars{10} --bs-trees 200 --threads $threads
raxml-ng --bootstrap --msa ${file_prefix}.T2.raxml.rba --prefix ${file_prefix}.T4 --threads $threads --seed 2 --bs-trees 100


## RAxML add bootstrapping values on the best scoring tree
raxml-ng --support --tree ${file_prefix}.T3.raxml.bestTree --bs-trees ${file_prefix}.T4.raxml.bootstraps --prefix ${file_prefix}.T5 --threads $threads 


cat ${file_prefix}*bootstraps > ${file_prefix}.T4.raxml.merged.bootstraps
raxml-ng --support --tree ${file_prefix}.T3.raxml.bestTree --bs-trees ${file_prefix}.T4.raxml.merged.bootstraps --prefix ${file_prefix}.T5 --threads $threads 




