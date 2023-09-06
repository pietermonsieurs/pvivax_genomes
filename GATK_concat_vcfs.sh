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


## OPTION 1: query the gvcf_combinegvcfs file, loop over the different samples, using 
## a wild character. However, this will not sort the vcf according to the position
## so a more supervised approach might be better (see below)
src_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master/
out_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat

for chrom in ${chrom_list[@]}
do
    echo $chrom
    bcftools concat ${src_dir}/PvP01_${chrom}_v1.*.vcf.gz -Oz -o ${out_dir}/PvP01_${chrom}_v1.vcf.gz
done



## OPTION 2: 
chrom_region_file=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master_merged/chrom_list_regions.csv
vcf_master_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master/
vcf_out_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/


for chrom in ${chrom_list[@]}
do
    ## create concat command
    bcf_concat="bcftools concat "

    ## extract the different regions out of the chromosome_region file
    while IFS=, read -r chrom_dummy start end
    do
        echo "${chrom} -> ${start}:${end}"
        bcf_concat+="${vcf_master_dir}/PvP01_${chrom}_v1.${start}.${end}.vcf.gz " 
        # echo $bcf_concat
    done <<< "$(grep "^${chrom}" ${chrom_region_file})"

    bcf_concat+=" -Oz -o ${vcf_out_dir}/PvP01_${chrom}_v1.vcf.gz"

    $bcf_concat
done


## quick quality check to see whether the number of samples in 
## each of the vcf files is OK
for vcf_file in *.vcf.gz; do bcftools query -l $vcf_file | wc; done
