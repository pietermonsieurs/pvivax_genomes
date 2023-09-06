module load BioTools

vcf_dir=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/
core_genome_bed_file=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/Pv_coregenome.bed


## filtering on the core genome if only filtering on quality
## of snps and indels has been performed, and are still split
## up per chromosome
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
    vcf_file_in=${vcf_dir}/filter/PvP01_${chrom}_v1.filtered.vcf.gz
    vcf_file_out=${vcf_dir}/filter_core/PvP01_${chrom}_v1.filtered.core.vcf.gz

    # bcftools index ${vcf_file_in}
    bcftools view -R ${core_genome_bed_file} ${vcf_file_in} -Oz -o ${vcf_file_out}

done


## filtering on the core genome in the case the already a filtering
## on the minor allele frequency has been performed

mafs=(0.005 0.010 0.050 0.100 0.200)
# mafs=(0.200)

for maf in ${mafs[@]}
do
    #maf=0.200
    core_genome_bed_file=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/Pv_coregenome.bed
    vcf_file_in=PvP01_v1.filtered.maf_${maf}.snps.vcf.gz
    vcf_file_out=PvP01_v1.filtered.maf_${maf}.snps.core.vcf.gz

    bcftools index ${vcf_file_in}
    bcftools view -R ${core_genome_bed_file} ${vcf_file_in} -Oz -o ${vcf_file_out}
done