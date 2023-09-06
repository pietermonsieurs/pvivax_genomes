#!/usr/bin/bash 

#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --job-name=QC_cov

module load BioTools
module load atools

## load all bam_files from the bam_master directory
source <(aenv --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/bam_master/bam_files.csv)

## minimum depth required. This is set to 5 as this is 
## also the value used by MalariaGen 
min_depth=5


## genome size: this is defined as the sum of all 
## core chromosomes, so partial contigs are removed, 
## and this number is hardcoded
genome_size=24250245

## bam_file input should be read from the input file
# bam_file=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/bam_master/SRR828528.bam
freq_coverage_file=${bam_file/.bam/.coverage}
freq_coverage_file=${freq_coverage_file/bam_master/bam_master_depth}
depth_file=${bam_file/.bam/.depth}
depth_file=${depth_file/bam_master/bam_master_depth}


## run the samtools depth command
samtools depth -a ${bam_file} > ${depth_file}

## calculate fraction of genome covered at least x times


# fraction=$(bc -l <<< "grep "^PvP01_" ${depth_file} | \
#      cut -f 3 -  | \
#      awk '$1>=$min_depth{c++} END{print c+0}' - )/${genome_size}")
fraction=$(grep "^PvP01_" $depth_file | awk -v x=${min_depth} '{if($3 >= x) covered++} END {total=24250245; print covered/total}'
)
mean=$(less ${depth_file} | awk '{sum+=$3} END {print sum/NR}')
median=$(less ${depth_file} | awk '{print $3}' | sort -n | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }')   
min=$(less ${depth_file} | awk 'BEGIN {min=1000000; max=0} {if ($3<min) min=$3; if($3>max) max=$3} END {print min}')
max=$(less ${depth_file} | awk 'BEGIN {min=1000000; max=0} {if ($3<min) min=$3; if($3>max) max=$3} END {print max}')



# fraction=$(bc -l <<< "$(samtools depth ${bam_file} | \
#     grep "^PvP01_" | \
#     cut -f 3 -  | \
#     awk '$1>=$min_depth{c++} END{print c+0}' - )/${genome_size}")

# samtools depth -a your_file.bam | awk -v x=10 '$3 >= x {total++; covered++} END {print "Fraction covered by at least", x, "reads:", covered/total}'
# samtools depth -a your_file.bam | awk -v x=10 '{if($3 >= x) covered++} {total++} END {print "Fraction covered by at least", x, "reads:", covered/total}'


# less $depth_file | awk -v x=5 '$3 >= x {total++; covered++} END {print "Fraction covered by at least", x, "reads:", covered, " - ", total, " - ", covered/total}'
# grep "^PvP01_" $depth_file | awk -v x=${min_depth} '{if($3 >= x) covered++} END {total=24250245; print covered/total}'

# mean=$(samtools depth -a your_file.bam | awk '{sum+=$3} END {print sum/NR}')
# median=$(samtools depth -a ${bam_file} | awk '{print $3}' | sort -n | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }')   

## extra sample name out of the bam_file
sample=${freq_coverage_file##*/}
sample=${sample/.coverage/}

## write down the output
echo "${sample},${fraction},${median},${mean},${min},${max}" > $freq_coverage_file

## clean up 
rm ${depth_file}

## run for all samples using atools
# cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/bam_master
# ls $PWD/*.bam > bam_files.csv
# vi /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/bam_master/bam_files.csv  ## add header bam_file
# module load atools
# sbatch --array $(arange --no_sniffer --data /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/bam_master/bam_files.csv) /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/bin/BWA_QC_percentage_mapped.sh