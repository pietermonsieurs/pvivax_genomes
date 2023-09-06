# Plasmodium vivax genomes

## Input data
* combined the fastq-file which are coming from multiple accession numbers: 
    * [inputdata_download_concate_multiple.py](inputdata_download_concate_multiple.py)
    * this script also moves some files between the source directories and some newly created. 
* overview of all sequencing data is given in the pvgenomes_master.xlsx excel file

## BWA and GATK
* Run BWA [BWA_run.sh](BWA_run.sh)
    * first map versus the human genome, and only take those reads not properly mapping
    * map remaining reads versus P. vivax P01 genome (version 46 of PlasmoDB)
* Remove samples with less than 50% of the genome covered by at least 5x reads
    * calculate percentage [BWA_QC_percentage_mapped.sh](BWA_QC_percentage_mapped.sh)
    * based on the percentages calculated, do filtering to remove low-quality genomes [BWA_QC_percentage_mapped_filter.sh](BWA_QC_percentage_mapped_filter.sh)


## location information
* Latitude and longitude information can be useful to finetune the geographical spread: python script for mapping: [inputdata_get_geolocation.py](inputdata_get_geolocation.py)
* map information on geographical map: show the total number of genomes available per country on a world wide map
    * update the R code [geomap_genomes.R](geomap_genomes.R)
    * now also plot circles based on coordinate (lat / long) give from metadata --> more finegrained view on the P vivax genomes in South-America
* same for the genomes in South-America, but now not one data point per country, but count given based on the availability of the longitude and latitude of the sampling. 
    
## Create vcf file for South-America only 
* create group with only the samples from the Americas
    * only select the samples coming from America (SAM) bases on sample IDs
    * [GATK_filter_bcftools_percontinent.sh](GATK_filter_bcftools_percontinent.sh)
* afterwards, concatenate the vcf files coming from different chromosomes. Use the code as already implmented in [GATK_filter_bcftools_concat.sh](GATK_filter_bcftools_concat.sh) where new code is added to run per continent



## Submission of the data to EVA repository of EBI
* chromosome names in the initial vcf file are taken from the PlasmoDB v46 genome. However, EBI requires mapping versus the genome assembly as available in RefSeq database: [EVA_update_chromosome_names.sh](EVA_update_chromosome_names.sh)
* concatenate all chromosomes into one overall vcf file: [EVA_bcftools_concat_chromosomes.sh](EVA_bcftools_concat_chromosomes.sh)
* remove the PL field in the vcf-fields as not requested by EVA: [EVA_bcftools_remove_PLfield.sh](EVA_bcftools_remove_PLfield.sh)
* EVA requires the BioSample accession number (not the run accession number), so run accession numbers need to be converted to BioSample accession number: [EVA_convert_runaccession_to_biosample.py](EVA_convert_runaccession_to_biosample.py)
* check validitiy of the vcf file using the vcf_validator tool of EBI: [EVA_vcf_validator.slurm](EVA_vcf_validator.slurm)




### procedure to add new data sets
* extract all meta data information from SRA using the ffq implementation
    * output a JSON file containing all experiments and runs related to that accession number [inputdata_download_from_SRA_ffq.sh](inputdata_download_from_SRA_ffq.sh)
    * afterwords, parse the fastq file URLs out of the JSON files: [inputdata_download_from_SRA_ffq.py](inputdata_download_from_SRA_ffq.py): this produces two files: one with additional information, one with a list of fastq files.
* list of fastq-files can be downloaded using: [inputdata_download_from_SRA.slurm](inputdata_download_from_SRA.slurm)
    * for sequencing data of Hugo Valdivia: additional script to concatenate fastq files from different runs [inputdata_prepareinput_valdivia.py](inputdata_prepareinput_valdivia.py). 
    * also added the FastQ file of sample of Erin to this directory (/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/fastq_valdivia) -> Pv21-19




## Analysis

### BWA & GATK
* BWA
    * first against human, and non-proper paired reads realigned to Pvivax genomes
    * [BWA_run.sh](BWA_run.sh) (including filtering of duplicated reads)
* GATK: 
    * run GATK directly per chrom, and write to directory per crhom in output diretory /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master/
    * script: [GATK_run_perchrom.slurm](GATK_run_perchrom.slurm)
    * afterwards do QC on the gvcf files: [GATK_run_perchrom_checkintegrity.sh](GATK_run_perchrom_checkintegrity.sh)
* combining gvcf files and do genotyping
    * !!OLD!! genotyping the combined gvcfs - !!!OLD!!! approach: 
        * combine all gvcf files: [GATK_combinegvcfs.slurm](GATK_combinegvcfs.slurm)
        * running the combined genotyping and convert the g.vcf file to a normal .vcf file [[GATK_genotype_combinegvcfs.slurm]([GATK_genotype_combinegvcfs.slurm] )
    * uses GenomicsDBImport to import the different gvcf files, and split into different chromosomes and regions
        * [GATK_genotype_genomicsDB_regions.py](GATK_genotype_genomicsDB_regions.py): look up the size of the chromosomes and define the different bins that will be used for the GenomicsDBImport
        * [GATK_merge_genomicsDBimport_prepare.sh](GATK_merge_genomicsDBimport_prepare.sh): create input files for the genomicsDBimport command, which mainly includes createing the mapping files, to map between samle_name and vcf files. You can use this to overwrite the sample name too.
        * [GATK_merge_genomicsDBImport.slurm](GATK_merge_genomicsDBImport.slurm): import per chromosome and per region, and do the import of the genome. Risk of out of memory error, so now running with 14 threads (to limit number of processes on on node), and with batch size of only 10 (instead of 50)
        * [GATK_genotype_genomicsDB.slurm](GATK_genotype_genomicsDB.slurm): genotype the databases imported via GenomicsDBImport
        * concatenate all vcf files per chromosome: [GATK_concat_vcfs.sh](GATK_concat_vcfs.sh)

* filtering of the SNP and indels files in different steps
    * apply the filtering procedures as suggested by GATK for indels and SNPs separately. Afterwards, merge all the data again: [GATK_filter_vcfs.sh](GATK_filter_vcfs.sh)
    * only select the core genome from the overall vcf file: [GATK_filter_bcftools_coregenome.sh](GATK_filter_bcftools_coregenome.sh)
    * filter out those samples with coverage lower than 50% covered at least 5x 
        * [BWA_QC_percentage_mapped.sh](BWA_QC_percentage_mapped.sh): creates one output file per genome, which should be concatenated. Afterwards, add this information to the meta data master file, and select those samples with a percentage smaller than 50%
        * [BWA_QC_percentage_mapped_filter.sh](BWA_QC_percentage_mapped_filter.sh): filter out those samples that do not have the minimum coverage. this is a copy-paste from the filtering in the excel file (locally)
            * includes also the option / bash command to remove those SNPs that are not genotype in the majority of the samples, for example exclude all positions with more than 10% missing genotypes. 
    * Old scripts but recycled in the new workflow
        * filter based on allele freqencies using bcftools: [GATK_filter_bcftools.sh](GATK_filter_bcftools.sh) --> different minor allele frequencies will be run, and applied on each chromosome file separately. 
        * concatenate the different chromosomes into one big file: [GATK_filter_bcftools_concat.sh](GATK_filter_bcftools_concat.sh) 


### Statistics
* Mapping of input data on geographical map
    * export from the master excel file
    * perform mapping using [geomap_genomes.R](geomap_genomes.R)
* get sequencing depth statistics based on the gvcf files of the different samples
    * [inputdata_getdepth.slurm](inputdata_getdepth.slurm)
    * extract the stats from the : [inputdata_getdepth_analysis.py](inputdata_getdepth_analysis.py)
    * plotting the depth stats: [inputdata_getdepth_plotting.R](inputdata_getdepth_plotting.R)


### Visualisation in 2D plot
* all data have been visualised in 2D plots with PCA, PCoA and UAMP in the Jupyter Notebook [UMAP.ipynb](UMAP.ipynb)
    * input file: PvP01_v1.filtered.maf_0.200.snps.core.vcf.gz
        * input contained SNPs with minimum allele frequency of 0.20. Can later on be decreased in RAM allows. 
        * also only core genome selected. Core genome file was downloaded from the MalariaGen website, and the core was saved as .bed file in the data directory on the supercomputer
* update: new analysis using PCA with PLINK
    * PCA analysis: [plink_ldpruning_and_file_conversion.sh](plink_ldpruning_and_file_conversion.sh)
    * this analysis also includes the LD pruning step
    * Plotting happens in same script as the admixture plotting script: [admixture_plotting.R](admixture_plotting.R)
    
### Admixture analysis
* admxiture software
    * https://speciationgenomics.github.io/ADMIXTURE/: requires installation of Admixture package
    * should be run on the output of the PLINK software, as this software puts the data in the correct input format
    * to be run for different number of K (= different populations) using cross-validation, and then select the one with the lowest error
    * [plink_admixture_CV.sh](plink_admixture_CV.sh): takes long time to run for larger K-values (up to 71 hours)
        * concatenate the different results by parsing out one line, containing CV error `grep 'CV error' *.log > admixture.Kvalues.CrossValidation.csv`
    * https://evomics.org/wp-content/uploads/2020/01/HowNotOverinterpretTutorialCesky2020.pdf
    
* PCA analysis using PLINK
    * crucial script: [plink_ldpruning_and_file_conversion.sh](plink_ldpruning_and_file_conversion.sh)
    * Should be prevented that there are samples with too many missing genotypes, and prevent that there are too many variants with missing genotypes in too many samples
        * removing variants which result frequently in missing genotype
            * either use the output files filter on the minor allele frequency with bcftools (0.005, 0.010, 0.050, 0.100, 0.200)
            * or use the --maf option in the PLINK2 software, eventually combined with the --geno flag where the percentage of missing genotypes per variant can be filtered out. 
        * removing samples with too many missing genotype: use the --mind flag in PLINK. For example, if set to 0.50 (so filter out samples with more than 50% missing genotypes), there is one sample that is removed (ERR2309687)
        * afterwards, remove the sample(s) from the vcf file using BCFtools: bcftoolf view -s ^ERR2309687 ....
* Phylogenetic tree 


IBD --> IsoRelate or PLINK


## TO DO
* update the download process. This could be optimised by using the ffq software (https://github.com/pachterlab/ffq) which allows to get the ftp-donwload link together with the file size and MD5 value, whic makes it much easier to follow up the download process
* get access to the correct location of the different latin american samples -> latitude and longitude
* integrate the additional 170 Peruvian samples (publication 2022)
* check 50% coverage of 5x for all genomes
* admixture analysis + finestructure/ChromPainter + ... 






        
