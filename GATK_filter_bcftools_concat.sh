module load BioTools

src_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/
cd $src_dir


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

## define all minor allele frequencies used for filtering
mafs=(0.005 0.010 0.050 0.100 0.200)
# mafs=(0.200)


## concat for the full list of variants (including indels)
for maf in ${mafs[@]}
do
    ## move to the working directory
    working_dir=${src_dir}/filter_core_coveragefilter_gtmissing_maf${maf}/
    cd ${working_dir}

    concat_command="bcftools concat -Oz "
    for chrom in ${chrom_list[@]}
    do
        concat_command+="PvP01_${chrom}_v1.filtered.core.coverage_filter_gtmissing.maf${maf}.vcf.gz "
        #echo $chrom
    done
    concat_command+="-o PvP01_v1.filtered.core.coverage_filter_gtmissing.maf${maf}.vcf.gz "
    echo $concat_command
    $concat_command
done


## concat for only SNPs
for maf in ${mafs[@]}
do
    ## move to the working directory
    working_dir=${src_dir}/filter_core_coveragefilter_gtmissing_maf${maf}/
    cd ${working_dir}

    concat_command="bcftools concat -Oz "
    for chrom in ${chrom_list[@]}
    do
        concat_command+="PvP01_${chrom}_v1.filtered.core.coverage_filter_gtmissing.maf${maf}.snps.vcf.gz "
        #echo $chrom
    done
    concat_command+="-o PvP01_v1.filtered.core.coverage_filter_gtmissing.maf${maf}.snps.vcf.gz "
    echo $concat_command
    $concat_command
done


#### do concat without filtering on minor allele frequency? 
cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/filter_core_coveragefilter_gtmissing/
concat_command="bcftools concat -Oz "
for chrom in ${chrom_list[@]}
do
    concat_command+="PvP01_${chrom}_v1.filtered.core.coverage_filter_gtmissing.vcf.gz "
    #echo $chrom
done
concat_command+="-o PvP01_v1.filtered.core.coverage_filter_gtmissing.vcf.gz "
echo $concat_command
$concat_command



##### PER CONTINENT #######

## concatenate the per-continent vcf files
cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_combinegvcf/subsets/

## first create a chromosome list by looping over the different 
## chromosomes
chrom_list=(01)
for chrom in {02..14}
do
    echo $chrom
    chrom_list+=(${chrom})
done

## all variants
concat_command="bcftools concat -Oz "
for chrom in ${chrom_list[@]}
do
    concat_command+="PvP01_${chrom}_v1.filtered.SAM.vcf.gz "
done
concat_command+="-o PvP01_v1.filtered.SAM.vcf.gz "
$concat_command

## all SNPs
concat_command="bcftools concat -Oz "
for chrom in ${chrom_list[@]}
do
    concat_command+="PvP01_${chrom}_v1.filtered.SAM.snps.vcf.gz "
done
concat_command+="-o PvP01_v1.filtered.SAM.snps.vcf.gz "
$concat_command


## all variants -- core
concat_command="bcftools concat -Oz "
for chrom in ${chrom_list[@]}
do
    concat_command+="PvP01_${chrom}_v1.filtered.SAM.core.vcf.gz "
done
concat_command+="-o PvP01_v1.filtered.SAM.core.vcf.gz "
$concat_command

## all SNPs -- core
concat_command="bcftools concat -Oz "
for chrom in ${chrom_list[@]}
do
    concat_command+="PvP01_${chrom}_v1.filtered.SAM.snps.core.vcf.gz "
done
concat_command+="-o PvP01_v1.filtered.SAM.snps.core.vcf.gz "
$concat_command