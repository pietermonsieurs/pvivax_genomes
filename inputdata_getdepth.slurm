#!/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --mem-per-cpu=4g
#SBATCH --time=0:30:00
#SBATCH --job-name=depth_vcf

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=pmonsieurs@itg.be


## load modules
module load atools
module load BioTools

## read in using atools the gvcf files
source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master/gvcf_files.csv)

## parameters
# gvcf_file=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master/FCH5YG3CCX2_L5_WHRDVIVqyuRAAHAA-92_NoHuman_VariantCalled_CoordSorted.g.vcf.gz
# gvcf_file=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master/19111936.g.vcf.gz
depth_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/depth_from_gvcf/


## create the filename for the depth file
gvcf_file=${gvcf_file}
file_prefix_full=${gvcf_file%.g.vcf.gz}
sample=${file_prefix_full##*/}
depth_file=${depth_dir}/${sample}.depth.csv

echo $gvcf_file
echo $depth_file

## run BCFtools
bcftools query -f '%CHROM,%POS[,%DP]\n' ${gvcf_file} > ${depth_file}

## run analysis directly on the depth_file using a python script
/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/inputdata_getdepth_analysis.py ${depth_file}

## run atools
# cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master
# ls $PWD/*.g.vcf.gz > gvcf_files.csv
# vi gvcf_files.csv ## add header gvcf_file
# module load atools
# sbatch --array $(arange --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master/gvcf_files.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/inputdata_getdepth.slurm
