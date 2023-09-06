module load BioTools

cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/filter_noPL/

## concatenatet the different vcf-files
bcftools concat -Ob PvP01_01_v1.filtered.no_PL.vcf.gz PvP01_02_v1.filtered.no_PL.vcf.gz PvP01_03_v1.filtered.no_PL.vcf.gz PvP01_04_v1.filtered.no_PL.vcf.gz PvP01_05_v1.filtered.no_PL.vcf.gz PvP01_06_v1.filtered.no_PL.vcf.gz PvP01_07_v1.filtered.no_PL.vcf.gz PvP01_08_v1.filtered.no_PL.vcf.gz PvP01_09_v1.filtered.no_PL.vcf.gz PvP01_10_v1.filtered.no_PL.vcf.gz PvP01_11_v1.filtered.no_PL.vcf.gz PvP01_12_v1.filtered.no_PL.vcf.gz PvP01_13_v1.filtered.no_PL.vcf.gz PvP01_14_v1.filtered.no_PL.vcf.gz PvP01_API_v1.filtered.no_PL.vcf.gz PvP01_MIT_v1.filtered.no_PL.vcf.gz -o PvP01_v1.filtered.no_PL.vcf.gz


## check the validity of the concatenated vcf file
/user/antwerpen/205/vsc20587/data/software/vcf-validator/vcf_validator_linux --input PvP01_v1.filtered.no_PL.vcf.gz --report text




## do the same but now for the correct number of genomes
cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/filter_core_coveragefilter_noPL

bcftools concat -Ob PvP01_01_v1.filtered.core.coverage_filter.no_PL.vcf.gz PvP01_02_v1.filtered.core.coverage_filter.no_PL.vcf.gz PvP01_03_v1.filtered.core.coverage_filter.no_PL.vcf.gz PvP01_04_v1.filtered.core.coverage_filter.no_PL.vcf.gz PvP01_05_v1.filtered.core.coverage_filter.no_PL.vcf.gz PvP01_06_v1.filtered.core.coverage_filter.no_PL.vcf.gz PvP01_07_v1.filtered.core.coverage_filter.no_PL.vcf.gz PvP01_08_v1.filtered.core.coverage_filter.no_PL.vcf.gz PvP01_09_v1.filtered.core.coverage_filter.no_PL.vcf.gz PvP01_10_v1.filtered.core.coverage_filter.no_PL.vcf.gz PvP01_11_v1.filtered.core.coverage_filter.no_PL.vcf.gz PvP01_12_v1.filtered.core.coverage_filter.no_PL.vcf.gz PvP01_13_v1.filtered.core.coverage_filter.no_PL.vcf.gz PvP01_14_v1.filtered.core.coverage_filter.no_PL.vcf.gz PvP01_API_v1.filtered.core.coverage_filter.no_PL.vcf.gz PvP01_MIT_v1.filtered.core.coverage_filter.no_PL.vcf.gz -o PvP01_v1.filtered.core.coverage_filter.no_PL.vcf.bgz


## check the validity of the concatenated vcf file --> better
## to run via a slurm script as this can take very long time
bcftools view -Ov PvP01_v1.filtered.core.coverage_filter.no_PL.vcf.gz -o PvP01_v1.filtered.core.coverage_filter.no_PL.vcf
/user/antwerpen/205/vsc20587/data/software/vcf-validator/vcf_validator_linux --input PvP01_v1.filtered.core.coverage_filter.no_PL.vcf --report text


## convert to a normal gzipped vcf-file
bcftools view -Oz -o PvP01_v1.filtered.core.coverage_filter.no_PL.vcf.gz  PvP01_v1.filtered.core.coverage_filter.no_PL.vcf.bgz
bcftools index --tbi PvP01_v1.filtered.core.coverage_filter.no_PL.vcf.gz
