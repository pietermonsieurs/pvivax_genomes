module load BioTools

## create new directory, and save in that directory the 
## manually created file with the old and new chromosome
## names
cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/filter_core_coveragefilter_noPL/
mkdir chrom_rename
cd chrom_rename
vi chrom_names.txt

## rename the vcf-file
bcftools annotate --rename-chrs chrom_names.txt ../PvP01_v1.filtered.core.coverage_filter.no_PL.vcf.gz -Oz  -o ./PvP01_v1.filtered.core.coverage_filter.no_PL.chrom_renamed.vcf.gz


## update three locations in the vcf file where there is no
## correspondance between the reference sequence at NCBI and the 
## sequence at PlasmoDB:

# Line 131877: Chromosome LT635614.2, position 144623, reference allele 'TA'
# does not match the reference sequence, expected 'tN'
# Line 131878: Chromosome LT635614.2, position 144625, reference allele
# 'CAAGCAGAAGACGGCATACGAG' does not match the reference sequence, expected
# 'NNNNNNNNNNNNNNNNNNNNNG'
# Line 131879: Chromosome LT635614.2, position 144643, reference allele
# 'CGAGATCATCAAGT' does not match the reference sequence, expected
# 'NNNGATCATCAAGT'

mkdir /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/filter_core_coveragefilter_noPL/chrom_rename_edit_ref

cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/filter_core_coveragefilter_noPL/chrom_rename_edit_ref/
cp ../chrom_rename/* ./
gunzip PvP01_v1.filtered.core.coverage_filter.no_PL.chrom_renamed.vcf.gz
vi PvP01_v1.filtered.core.coverage_filter.no_PL.chrom_renamed.vcf ## edit locations - see above
wc PvP01_v1.filtered.core.coverage_filter.no_PL.chrom_renamed.vcf ## should be 1837232

gzip PvP01_v1.filtered.core.coverage_filter.no_PL.chrom_renamed.vcf

## you need a blocked gzip file (to be created with bcftools) to be
## able to run tabix on the file
module load BioTools
bcftools view -Ob PvP01_v1.filtered.core.coverage_filter.no_PL.chrom_renamed.vcf -o PvP01_v1.filtered.core.coverage_filter.no_PL.chrom_renamed.vcf.bgz


## create .tbi file
module load BioTools
tabix PvP01_v1.filtered.core.coverage_filter.no_PL.chrom_renamed.vcf.bgz


