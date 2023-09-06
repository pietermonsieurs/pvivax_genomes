module load BioTools

samples_to_filter=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/bam_master_depth/coverage_5x.too_low.csv
vcf_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/

## create chromosome list
chrom_list=(01)
for chrom in {02..14}
do
    echo $chrom
    chrom_list+=(${chrom})
done

chrom_list+=("API")
chrom_list+=("MIT")
echo ${chrom_list[@]}

for chrom in ${chrom_list[@]}
do
    echo $chrom
    vcf_file_in=${vcf_dir}/filter_core/PvP01_${chrom}_v1.filtered.core.vcf.gz
    vcf_file_out=${vcf_dir}/filter_core_coveragefilter/PvP01_${chrom}_v1.filtered.core.coverage_filter.vcf.gz
    
    bcftools index ${vcf_file_in}
    bcftools view -S ^${samples_to_filter} ${vcf_file_in} -Oz -o ${vcf_file_out}

done


for chrom in ${chrom_list[@]}
do
    echo $chrom
    vcf_file_in=${vcf_dir}/filter_core_coveragefilter/PvP01_${chrom}_v1.filtered.core.coverage_filter.vcf.gz
    vcf_file_out=${vcf_dir}/filter_core_coveragefilter_gtmissing/PvP01_${chrom}_v1.filtered.core.coverage_filter_gtmissing.vcf.gz
    
    bcftools index ${vcf_file_in}
    bcftools view ${vcf_file_in} \
        -i 'F_MISSING<0.10' \
        -Oz -o ${vcf_file_out}
done



## debug / QC vcf -files
for chrom in ${chrom_list[@]}
do
    echo $chrom
    vcf_file_in=${vcf_dir}/filter/PvP01_${chrom}_v1.filtered.vcf.gz
    bcftools query -l ${vcf_file_in} | wc
done

cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master
for chrom in ${chrom_list[@]}
do
    echo $chrom
    for vcf_file in PvP01_${chrom}_v1*.vcf.gz
    do
        bcftools query -l ${vcf_file} | wc
    done
done