module load BioTools
module load Java
module load GATK



## first create a chromosome list by looping over the different 
## chromosomes
chrom_list=(01)
# for chrom in {02..14}
# chrom_list=(07)
for chrom in {02..14}
do
    echo $chrom
    chrom_list+=(${chrom})
done

chrom_list+=("API")
chrom_list+=("MIT")
echo ${chrom_list[@]}


## check integrity of the g.vcf file
# src_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_perchrom/
src_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master/
meta_data_file=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/PvGenomes_master.txt
bam_master_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/bam_master/



## select column 21 (search_name_print) and not the search_name
## column as this still contains slashes (/). Remove the first elelement
## from the header as this is the column name
sampleIDs=($(cut -f 21 $meta_data_file))
# unset sampleIDs[0]

## loop over the samples extracted from the master file and check
## whether the g.vcf file exists. If the same for all chromosomes, only
## run for one chrom (e.g. chrom01) and select the sample ID to rerun GATK
for sampleID in "${sampleIDs[@]}"
do
    # echo $sampleID
    for chrom in ${chrom_list[@]}
    do
        # echo $chrom
        chrom_dir=${src_dir}/PvP01_${chrom}_v1/

        gvcf_file=${chrom_dir}/${sampleID}.PvP01_${chrom}_v1.g.vcf.gz

        # gatk ValidateVariants -V $gvcf_file
        # if ! gzip -t $gvcf_file
        if ! test -f "$gvcf_file"
        # then echo -e "${sampleID}\tfile ${gvcf_file} is NOT ok"
        then echo ${bam_master_dir}/${sampleID}.bam
        # else echo "file ${gvcf_file} ok"
        fi
    done
done > debug.out


## extract the problematic samples from the debug file, as 
## the sample ID is in the first column
samples_missed=($(cut -f 1 debug.out))
bam_master_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/bam_master/
for sampleID in "${samples_missed[@]}"
do 
    bam_file=${bam_master_dir}/{$sampleID}.bam
    echo $bam_file
done

