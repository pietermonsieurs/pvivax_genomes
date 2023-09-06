## vcf-validator can direct be downloaded from the github 
## page of EBI
cd /user/antwerpen/205/vsc20587/data/software/vcf-validator
wget https://github.com/EBIvariation/vcf-validator/releases/download/v0.9.4/vcf_validator_linux

## run for some of the chromosomes
./vcf_validator_linux --input  ~/scratch/plasmodium_pvgenomes/results/vcf_master_concat/PvP01_01_v1.vcf.gz

## exporting as text gives more information on the error message. Problem
## lies with the PL field that is not always reporting the correct number of 
## fields. Removing this from the vcf file with bcftools solves the problem. 
## can be donw with EVA_bcftools_remove_PLfield.sh
./vcf_validator_linux --input  ~/scratch/plasmodium_pvgenomes/results/vcf_master_concat/PvP01_API_v1.no_PL.vcf --report text



## debug message --> write separate python script
