# Plasmodium vivax genomes

## Overview of Data

### Daaron et al. paper 2021
Overview of the different genomes is given in the supplementary table: https://www.science.org/doi/suppl/10.1126/sciadv.abc3713/suppl_file/abc3713_table_s1.xlsx --> now saved in /Users/pmonsieurs/programming/plasmodium_pvgenomes/data/2016.Pearson.list_of_accesionnumbers.PvGV_May2016_sample_data_2.xlsx
* This database contains data from 12 different sequencing projects, resulting in 1154 P. vivax isolates. The 12 project numbers are PRJEB10888, PRJEB2140, PRJNA175266, PRJNA240356 , PRJNA240452, PRJNA240531, PRJNA271480 , PRJNA284437 , PRJNA350554 , PRJNA420510 , PRJNA432819 , PRJNA65119 
* Within the article itself, 20 additional genomes from Africa (and isolated from human) are sequences: Mauritania, n = 14; Ethiopia, n = 3; and Sudan, n = 3. Those seem to be grouped in project number PRJNA720520. This project number also contains 10 genome sequences obtained from apes and mosquitoes, which should be omitted. 
* in total: 1154 + 20 new in human + 10 in other species = 1184 genomes <--> in the excel of the paper only 1181 genomes, so some of them are missing. 
* list of Daaron et al. also contains Pvivax like genomes. They should be removed from the list, as they are not infecting human hosts. Those are 27 WGS sequences (last lines in the excel file). 
* some of the accession numbers in Daaron et al are not the same as in the data set of Juie, because Daaron et al sometimes splits them up in different samples, while they should be combined into one sample
    * Julie has made this distinction. summarized now in "concatenated_list" in /Users/pmonsieurs/programming/plasmodium_pvgenomes/data/julie/Metadata.xlsx

* make comparison between two databases (Daaron vs Julie): [inputdata_compare_db.py](inputdata_compare_db.py)
    * start from exel file from Daaron et al and compare to the excel file of Julie
    * two different approach: first check accession numbers in both "main" work sheets, and - if no hit found - check the accession number with the concatenated accession numbers of Julie (second work sheet)
    * most of the genomes added by Daaron, were already in the list of Juli. Some of them were missed by Julie, and are now added in a separate excel worksheet containing the accession numbers that still needs to be downloaded: Excel file in data directory: 2021.Daron.SupplTableS1.abc3713_table_s1.comparison_with_database_Julie.xlsx --> sheet: samples_to_add_to_Julie. Total number of samples to be added = 267 genomes (8 will be removed - see below)
        * missed publications:
            * Shen et al. 2018: China, n=7
            * Popovici et al. 2018: Cambodia n=55
            * Pearson et al. 2016: Cambodia n=11
            * Delgado-ratto et al. 2016: Brazil, n=28
        * new publications:
            *  Daaron et al 2021: Ethiopia n=3, Mauritania n=14, Sudan n=3
            *  Benavente et al 2021: Afghanistan n=22, Bangladesh n=1, Brazil n=27, Eritrea n=12, Ethiopia n=5, Guyana n=3, India n=34, Pakistan n=32, Sudan n=5, Philipinnes n=1, Uganda n=3
        * new in-house data or collab (!! not included yet in the 267 genomes list)
            * Erin: 15 samples from clinic (different countries, mainly Africa)
            * data from Marcello: n= ?? --> contact taken with Thais Crippa to transfer the data
    * some of the new data are single-end data? 
        * in total 267 genomes can be added, however only 259 have a *_1.fasta.gz extension which means 8 are missing.
        * single-end data + they end up to be IonTorrent data: SRR5099324.fastq.gz  SRR5099325.fastq.gz  SRR5099326.fastq.gz  SRR5278291.fastq.gz  SRR5278293.fastq.gz  SRR5278295.fastq.gz  SRR5278298.fastq.gz  SRR5278301.fastq.gz. Remarkably, they ended up in the final .vcf file of Julie --> remove from master file

### MalariaGen paper (2022)
* based on file: /Users/pmonsieurs/programming/plasmodium_pvgenomes/data/MalariaGen_PV_samples.xlsx, which is a copy paste from the meta data from the paper. 
* Which samples added? Different steps in selection: 
	* Select only the samples where the counterpart in Julie is empty (column P)
	* Only select those samples that survived the QC criteria, i.e. exclusion reason = Analysis_set
    * Only select samples for which the ENA accession number is given.
    * removal of the samples of Huppalo, as those 6 samples are already in Julie's dataset
    * 535 genomes are retained
* all those samples are summarised in the last tab of the excel file
    * tab: genomes_add_to_master
    * in total 571 datasets can be downloaded (should be 575)
        * for 20 samples [row 101 + row 454-472], three accession numbers are found = 40 additional dataset
        * 571  - 40 = 531 datasets left <--> should be 535: 4 datasets are missed between row 101 and 200
    * put accession numbers into SRA-explorer, in batches of ~ 100, and put in "saved dataset"
        * when creating download links: 567 download links found, i.e. another 4 accession numbers lost
* download fastq files based on the download links created with SRA explorer
    * !! OBSOLETE!! trying to download the data by using atools combined with wget [inputdata_download_from_SRA.slurm](inputdata_download_from_SRA.slurm), however this leads to frequent crasehs. Seems that only one connection can be made between a node and the EBI server, all other connections are refused.
    * using a basic solution: splitting up the big URL file into different smaller files [input_download_from_SRA_batches.sh](input_download_from_SRA_batches.sh)
        * run the download using a simple for loop for each of the different batches
* track the missing accession numbers: some of the accession numbers are not found when running SRA explorer, and for others the URLs are not found
    * find missing accession numbers: [inputdata_download_find_missing_malariagen.py](inputdata_download_find_missing_malariagen.py)
    * after running again this list of 8 accession numbers (ERR2299660,ERR925422,ERR925423,ERR925429,ERR1035495,ERR2309737,ERR2299732,ERR2299718), 7 are recovered but one is still missing: ERR2299660
    * after downloading all the ftp-URLs, total number of files downloaded is 1148 = 574 datasets of paired fastq file
        * 535 genomes of which 20 with 3 accession numbers = 575 R1 fastq files
        * only 574 R1 fastq files found as one accession number (ERR2299660 - sampleID PY0130-C) is missing. FastQ is also not available when chekcing online: https://www.ebi.ac.uk/ena/browser/view/ERR2299660 --> fastq column: "unavailable"
* combined the fastq-file which are coming from multiple accession numbers: 
    * [inputdata_download_concate_multiple.py](inputdata_download_concate_multiple.py)
    * this script also moves some files between the source directories and some newly created. 
    * after moving concatenated fastqfiles again to the source directory, and afterwards check again whether then number of R1 reads is sufficient --> 535 genomes with one accession number continously failing = 534 genomes



### Map metadata julie to VCF files ###
* The Metadata.xlsx file of Julie contains all metadata information (country, accessionnumbers(s), etc.) from the different samples. However, there is not always a direct link between the metadata file and the name of the vcf-file. Moreover, in the overal / combined vcf file, the name of the vcf file has been replaced by a more meaningful name e.g. Pearson_Vietnam_3
    * in most cases, the accesssion number is the link between both files
    * however, when multiple accession numbers are linking to the same sample, sometimes a new name is given to the .vcf file 
    * also, for the inhouse samples the naming cannot follow the accession number rule, as they have not yet been uploaded to NCBI-SRA
* once the metadata file has been updated, a newly compiled metadata can be produced. 
* !! data of Pearson et al: old accession numbers, not used anymore. Difficult to map to new accession numbers, which was used in subsequent steps by Julie but no mapping. Paper of G. Spath on GIP brought the solution: https://www.biorxiv.org/content/10.1101/2021.06.15.448580v1
* Metadata.xlsx adapted: added additional columsn to store the mapping to the .vcf and the new sample names used in th merged .vcf file
    * used MetaData.xslx as source for a new metadata file



### create new master file
* new file created: PvGenomes_Master.xlsx
    * remove all samples with NA value in the vcf_old (this means those samples do not have a vcf file in the output of Julie)
        * 6 NAs deleted from sample data as no counterpart in vcf files
    * still one more in sample file than in vcf files --> duplicate? 

### new data
* New data recovered from the Daaron paper + new papers published in 2021 / end of 2020
    * All the fastq-file R1 files are used as input for running BWA. BWA is run in 2 times, one time versus the human genome, and the second time using the Pvivax reference genome using all reads that did not map to the human genome. For now, those data have been stored in the directory 
    * [BWA_run.sh](BWA_run.sh): run BWA based on the first read fastq file
    * Some of the input data are too big to be run within one hour, and should be moved to 24h queue. This initially causes a crash of the run, and recovering of the missed accession numbers using [BWA_recover_crashed.py](BWA_recover_crashed.py). Afterwards rerun the [BWA_run.sh](BWA_run.sh) but with a time of 23:59:00.
* New genomes Erin from clinic
    * make symbolic links based on the internal ID with the genomeScan IDs, and make sure to continue on the syntax of using _1.fastq.gz and _2.fastq.gz: [inputdata_prepareinput_erin.py](inputdata_prepareinput_erin.py)
    * once symbolic links are made, BWA can be run with the script developed before [BWA_run.sh](BWA_run.sh)
* New genomes from Brazil (Marcelo)
    * convert the unaligned bam-files to fastq files: [inputdata_prepareinput_marcelo_bam2fastq.sh](inputdata_prepareinput_marcelo_bam2fastq.sh)
    * make symbolic links [inputdata_prepareinput_marcelo_makelinks.sh](inputdata_prepareinput_marcelo_makelinks.sh)
* Data from Julie:
    * move all data from Julie (ordered in different subdirectories) into one flat directory: [inputdata_julie_filemanagement.sh](inputdata_julie_filemanagement.sh)
* new data from MalariaGen publication (2022)
    * [BWA_run.sh](BWA_run.sh): run BWA based on the first read fastq file


### reschedule files and MarkDuplicates
* MarkDuplicates has not been run for all samples. Moreover, the bam files of Julie are needed to do a QC of the different genomes (requiring that at least 50% of the genome is covered 5x)
    * map the different .bam files in the results directory to the correct sample: [inputdata_julie_map_metadata_to_bam.py](inputdata_julie_map_metadata_to_bam.py)
    * based on this output [/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/mapping_masterdata_to_bam.csv] move all bam files to a bam_master directory
        * for the files where no markDuplicates has been run: run mark duplicates and save output in this directory
        * [GATK_markdups.sh](GATK_markdups.sh)
        * afterwards, check whether all sampleIDs are represented in the bam_master directory [GATK_markdups_checkfiles.py](GATK_markdups_checkfiles.py) --> can be recycled for checking all sampleIDs (also the ones from Julie)
* for the .bam files of Julie received from Katlijn (where MarkDuplicates is alreayd performed), the files are moved and renamed to the bam_master directory
    * move bam files to bam_master: [inputdata_julie_movebamfiles.py](inputdata_julie_movebamfiles.py). Also copy the .bai file (index files)
    * QC check: when checking bai-files, some were missing. Trying to redo them, and running `samtools quickcheck -v *.bam` showed that some .bam files were empty (received from Katlijn)
    * for those 64 files, re-copy from Katlijn, and add to the master_bam directory: [inputdata_julie_movebamfiles_missing.py](inputdata_julie_movebamfiles_missing.py)
    * afterwards run sanity check: [GATK_markdups_checkfiles.py](GATK_markdups_checkfiles.py) -> contains 1533 bam-files, only Pv21-19 is missing
    * also run again `samtools quickcheck -v *.bam`
    

### location information
* Latitude and longitude information can be useful to finetune the geographical spread
    * try to move the information in the MalariaGen paper to the master excel-file
    * metadata of MalariaGen paper downloaded from https://www.malariagen.net/data/open-dataset-plasmodium-vivax-v4.0 into /Users/pmonsieurs/programming/plasmodium_pvgenomes/data/metadata_2022_malariagen.xlsx
    * mapping can be done based on the sample_ids
    * python script for mapping: [inputdata_get_geolocation.py](inputdata_get_geolocation.py)
    * merge data with MalariaGen, and use this information for longitude and latitude. 
        * if no match with MalariaGen, copy paste from MalariaGen, or look up on web
        * indicated with separate column whether inherited from MalariaGen or manually added
        * also: sometimes different positions merges, e.g. the different Peru samples
* map information on geographical map
    * update the R code [geomap_genomes.R](geomap_genomes.R)
    * now also plot circles based on coordinate (lat / long) give from metadata --> more finegrained view on the P vivax genomes in South-America
    
### Create subgroups 
* create group with only the samples from the Americas (+ maybe Americas + Africa?)
    * only select the samples coming from America (SAM) bases on sample IDs
    * [GATK_filter_bcftools_percontinent.sh](GATK_filter_bcftools_percontinent.sh)
* afterwards, concatenate the vcf files coming from different chromosomes. Use the code as already implmented in [GATK_filter_bcftools_concat.sh](GATK_filter_bcftools_concat.sh) where new code is added to run per continent


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
    * [BWA_run.sh](BWA_run.sh)
    * to be added: script to move bam-files from temporary directory to master directory!!
    * check the number of missing samples using [BWA_check_missing.py](BWA_check_missing.py)
        * 20230109: Pv21-19 missing (sample Peru from Katlijn)

* GATK: 
    * Old approach (obsolete!!): 
        * [GATK_run.slurm](GATK_run.slurm): create one big g.vcf file per sample
        * split one gvcf file over a per chromosome file, to make merging and genotyping of the different files easier
            * [GATK_split_gvcf.slurm](GATK_split_gvcf.slurm)
            * Results in one directory per chromosome, where each subsequent analysis is performed on a per chromosome-basis
    * New approach: 
        * remove duplidates from the bam file [GATK_markdups.sh](GATK_markdups.sh)
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






        
