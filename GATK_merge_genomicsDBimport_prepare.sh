## script to create map file for the GenomicsDB import

module load BioTools


## first create a chromosome list by looping over the different 
## chromosomes
chrom_list=(01)
for chrom in {02..14}
do
    echo $chrom
    chrom_list+=(${chrom})
done

chrom_list+=("API")
chrom_list+=("MIT")
echo ${chrom_list[@]}


## query the gvcf_master file, loop over the different samples
## also put in a list. Either you read the sample names from the 
## gvcf file, and use them again as sample names, or you use
## the mapping file to rename the samples (search_name_print column)
## First part of the code if for extracting the sample name from 
## the gvcf file, second part is for renaming samples. PART II is the 
## preferred way!!

src_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master/
merge_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master_merged/
meta_data_file=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/PvGenomes_master.txt


## PART 1: extract sample name from GVCF -- OBSOLETE
# src_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_perchrom/
    
for chrom in ${chrom_list[@]}
do
    echo $chrom
    chrom_dir=${src_dir}/PvP01_${chrom}_v1/
    cd $chrom_dir

    for gvcf_file in *.g.vcf.gz
    do 
        sample=$(bcftools query -l ${gvcf_file})
        echo -e "${sample}\t${gvcf_file}"
    done > ${merge_dir}/mapping_file.PvP01_${chrom}_v1.map
done



## PART 2: rename sample names using the name of the file
sampleIDs=($(cut -f 21 $meta_data_file))
unset sampleIDs[0]

## loop over the samples extracted from the master file and check
## whether the g.vcf file exists. If the same for all chromosomes, only
## run for one chrom (e.g. chrom01) and select the sample ID to rerun GATK
for chrom in ${chrom_list[@]}
do
    echo $chrom
    mapping_file=${merge_dir}/mapping_file.PvP01_${chrom}_v1.map
    for sampleID in "${sampleIDs[@]}"
    do
        gvcf_file=${src_dir}/PvP01_${chrom}_v1/${sampleID}.PvP01_${chrom}_v1.g.vcf.gz
        echo -e "${sampleID}\t${gvcf_file}"
    done > ${mapping_file}
done




## create a subset of the map files to check how faste
## the merging goes

subset=100
cd $merge_dir
for map_file in *_v1.map
do 
    echo $map_file
    map_out=${map_file/_v1.map/_v1.test.map}
    echo ${map_out}
    # head -n ${subset} ${map_file} >  ${map_out
done