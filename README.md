# Plasmodium vivax genomes

## Input data
* combined the fastq-file which are coming from multiple accession numbers: 
    * [inputdata_download_concate_multiple.py](inputdata_download_concate_multiple.py)
    * this script also moves some files between the source directories and some newly created. 
* overview of all sequencing data is given in the pvgenomes_master.xlsx excel file

 


### BWA & GATK
* Run BWA [BWA_run.sh](BWA_run.sh)
    * first map versus the human genome, and only take those reads not properly mapping
    * map remaining reads versus P. vivax P01 genome (version 46 of PlasmoDB)
    * Remove samples with less than 50% of the genome covered by at least 5x reads
        * calculate percentage [BWA_QC_percentage_mapped.sh](BWA_QC_percentage_mapped.sh)
        * based on the percentages calculated, do filtering to remove low-quality genomes [BWA_QC_percentage_mapped_filter.sh](BWA_QC_percentage_mapped_filter.sh)
* GATK: 
    * run GATK directly per chrom, and write to directory per crhom in output diretory /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master/
    * script: [GATK_run_perchrom.slurm](GATK_run_perchrom.slurm)
    * afterwards do QC on the gvcf files: [GATK_run_perchrom_checkintegrity.sh](GATK_run_perchrom_checkintegrity.sh)
* combining gvcf files and do genotyping: uses GenomicsDBImport to import the different gvcf files, and split into different chromosomes and regions
    * [GATK_genotype_genomicsDB_regions.py](GATK_genotype_genomicsDB_regions.py): look up the size of the chromosomes and define the different bins that will be used for the GenomicsDBImport
    * [GATK_merge_genomicsDBimport_prepare.sh](GATK_merge_genomicsDBimport_prepare.sh): create input files for the genomicsDBimport command, which mainly includes createing the mapping files, to map between samle_name and vcf files. You can use this to overwrite the sample name too.
    * [GATK_merge_genomicsDBImport.slurm](GATK_merge_genomicsDBImport.slurm): import per chromosome and per region, and do the import of the genome. Risk of out of memory error, so now running with 14 threads (to limit number of processes on on node), and with batch size of only 10 (instead of 50)
    * [GATK_genotype_genomicsDB.slurm](GATK_genotype_genomicsDB.slurm): genotype the databases imported via GenomicsDBImport
    * concatenate all vcf files per chromosome: [GATK_concat_vcfs.sh](GATK_concat_vcfs.sh)
* filtering of the SNP and indels files in different steps
    * apply the filtering procedures as suggested by GATK for indels and SNPs separately. Afterwards, merge all the data again: [GATK_filter_vcfs.sh](GATK_filter_vcfs.sh)
    * only select the core genome from the overall vcf file: [GATK_filter_bcftools_coregenome.sh](GATK_filter_bcftools_coregenome.sh)
    * Old scripts but recycled in the new workflow
        * filter based on allele freqencies using bcftools: [GATK_filter_bcftools.sh](GATK_filter_bcftools.sh) --> different minor allele frequencies can be run, and applied on each chromosome file separately. 
        * concatenate the different chromosomes into one big file: [GATK_filter_bcftools_concat.sh](GATK_filter_bcftools_concat.sh) 
* Create vcf file for South-America only 
    * create group with only the samples from the Americas
        * only select the samples coming from America (SAM) bases on sample IDs
        * [GATK_filter_bcftools_percontinent.sh](GATK_filter_bcftools_percontinent.sh)
    * afterwards, concatenate the vcf files coming from different chromosomes. Use the code as already implmented in [GATK_filter_bcftools_concat.sh](GATK_filter_bcftools_concat.sh) where new code is added to run per continent


## location information
* Latitude and longitude information can be useful to finetune the geographical spread: python script for mapping: [inputdata_get_geolocation.py](inputdata_get_geolocation.py)
* map information on geographical map: show the total number of genomes available per country on a world wide map
    * update the R code [geomap_genomes.R](geomap_genomes.R)
    * now also plot circles based on coordinate (lat / long) give from metadata --> more finegrained view on the P vivax genomes in South-America
* same for the genomes in South-America, but now not one data point per country, but count given based on the availability of the longitude and latitude of the sampling. 


## Population genomics
* analysis using PCA with PLINK
    * PCA analysis: [plink_ldpruning_and_file_conversion.sh](plink_ldpruning_and_file_conversion.sh)
    * this analysis also includes the LD pruning step
    * Plotting happens in same script as the admixture plotting script: [admixture_plotting.R](admixture_plotting.R)
* admixture software
    * [plink_admixture_CV.sh](plink_admixture_CV.sh): takes long time to run for larger K-values (up to 71 hours)
        * concatenate the different results by parsing out one line, containing CV error `grep 'CV error' *.log > admixture.Kvalues.CrossValidation.csv`
    * plot CV error to identify optimal number of K's: [plink_admixture_CV_plotting.R](plink_admixture_CV_plotting.R)
* phylogenetic trees
    * created with code [RAxML.sh](RAxML.sh)
    * visualisation using ggtree [RAxML_visualize_trees_ggtree.R](RAxML_visualize_trees_ggtree.R)
* calculation of the fws values using the moimix package: [fws_moimix.R](fws_moimix.R)
* visualisation of population genomics approach:
    * overview of the different genomes sampled worldwide: [geomap_genomes.R](geomap_genomes.R) (also subset for South America included)
    * admixture plotting and PCA plotting: [plink_admixture_plotting.R](plink_admixture_plotting.R)
    * plotting RAxML tree: [RAxML_visualize_trees_ggtree.R](RAxML_visualize_trees_ggtree.R)
    * bubble plot of admixture analysis: [geomap_admixure_bubbleplot_piechart.R](geomap_admixure_bubbleplot_piechart.R)
    * geographical map of South America with pie chart indicating the admixture populations: [geomap_admixture_piechart.R](geomap_admixture_piechart.R)


## Submission of the data to EVA repository of EBI
* chromosome names in the initial vcf file are taken from the PlasmoDB v46 genome. However, EBI requires mapping versus the genome assembly as available in RefSeq database: [EVA_update_chromosome_names.sh](EVA_update_chromosome_names.sh)
* concatenate all chromosomes into one overall vcf file: [EVA_bcftools_concat_chromosomes.sh](EVA_bcftools_concat_chromosomes.sh)
* remove the PL field in the vcf-fields as not requested by EVA: [EVA_bcftools_remove_PLfield.sh](EVA_bcftools_remove_PLfield.sh)
* EVA requires the BioSample accession number (not the run accession number), so run accession numbers need to be converted to BioSample accession number: [EVA_convert_runaccession_to_biosample.py](EVA_convert_runaccession_to_biosample.py)
* check validitiy of the vcf file using the vcf_validator tool of EBI: [EVA_vcf_validator.slurm](EVA_vcf_validator.slurm)


        
