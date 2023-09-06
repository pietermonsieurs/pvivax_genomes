
plink=/user/antwerpen/205/vsc20587/data/software/plink_1.90/plink
plink2=/user/antwerpen/205/vsc20587/data/software/plink_2.00/plink2
admixture=/user/antwerpen/205/vsc20587/data/software/admixture_linux-1.3.0/admixture


## Pv genomes WORLD 
maf=0.005
cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/plink/
vcf_file=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_combinegvcf/PvP01_v1.filtered.maf_${maf}.snps.core.vcf.gz
out_prefix=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/plink/PvP01_v1.filtered.maf_${maf}.snps.core

## Pv genomes SAM
cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/plink/
vcf_file=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/filter_core_coveragefilter_gtmissing_maf0.005/PvP01_v1.filtered.core.coverage_filter_gtmissing.maf0.005.snps.SAM.vcf.gz
out_prefix=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/plink/PvP01_v1.filtered.core.coverage_filter_gtmissing.maf0.005.snps.SAM

## Pv genomes SAM but without outliers 
vcf_file_outliers=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/plink/SAM/PvP01_v1.no_outliers.vcf.gz
bcftools view -s ^SM-IA1D5,SM-IA1D7,SM-IA1D8,Panama_458,SRP046126 $vcf_file -Oz -o ${vcf_file_outliers}
out_prefix=/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/plink/SAM/PvP01_v1.no_outliers
vcf_file=${vcf_file_outliers}
## PLINK -- version 1.9
# $plink --vcf $vcf_file --recode --allow-extra-chr --double-id --out test.ped 

# # https://www.cog-genomics.org/plink/2.0/data#set_all_var_ids

 
# $plink --file test.ped --allow-extra-chr --indep-pairwise 50 5 0.5
# $plink --file test.ped --allow-extra-chr --extract plink.prune.in --make-bed --out pruneddata
# $plink --bfile pruneddata --allow-extra-chr --write-snplist


## PLINK2
# https://core.ac.uk/download/pdf/286682225.pdf

## converting files to other file formats using the 
## pgen option
$plink2 --vcf $vcf_file \
    --allow-extra-chr \
    --out ${out_prefix} \
    --make-pgen \
    --threads 1 \
    --set-all-var-ids @:# # \
    #--maf 0.20


## do LD pruning: remove those SNPs that are too closely
## related to each other
# $plink2 binary_fileset --allow-extra-chr --indep-pairwise 50 5 0.5
$plink2 --indep-pairwise 50 5 0.5 \
    --pfile ${out_prefix} \
    --allow-extra-chr \
    --threads 1 


## filtering out the positions after LD pruning and
## create again a pgen-like file
$plink2 --pfile ${out_prefix} \
    --extract plink2.prune.in \
    --allow-extra-chr \
    --out ${out_prefix}.ldpruned \
    --threads 1 \
    --make-pgen


## run PCA on the LD pruned pgen-file
$plink2 --pfile ${out_prefix}.ldpruned \
    --pca biallelic-var-wts 5 \
    --out ${out_prefix}.pca_results \
    --allow-extra-chr \
    --threads 1 

## try to do LD Pruning but export to a bed file?
$plink2 --pfile ${out_prefix} \
    --extract plink2.prune.in \
    --allow-extra-chr \
    --out ${out_prefix}.ldpruned \
    --threads 1 \
    --make-bed


## option 2: convert the pgen file after removal of LD-pruned
## positions to a bed file using --make-bed
$plink2 --pfile ${out_prefix}  \
    --extract plink2.prune.in \
    --allow-extra-chr \
    --out ${out_prefix}.ldpruned.filtered \
    --threads 1 \
    --mind 0.50 \
    --make-bed

## not sure anymore why this needs to be done?
awk '{$1="0";print $0}' ${out_prefix}.ldpruned.filtered.bim > ${out_prefix}.ldpruned.filtered.bim.tmp
mv ${out_prefix}.ldpruned.filtered.bim.tmp ${out_prefix}.ldpruned.filtered.bim

$admixture ${out_prefix}.ldpruned.filtered.bed 2  > test.log
