# install.packages('remotes')
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("bahlolab/moimix", build_vignettes = TRUE)

## load the required libraries
library(moimix)
library(SeqArray)
library(openxlsx)
library(ggplot2)

## parameter settings. PLINK filtered means that the same samples
## are removed from the database as in the PLINK software
src_dir = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/gvcf_combinegvcf/'
src_name = 'PvP01_v1.filtered.maf_0.200.snps.core'
src_name = 'PvP01_v1.filtered.maf_0.200.snps.core.plink_filtered'
src_name = 'PvP01_v1.filtered.maf_0.100.snps.core'

src_dir = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/vcf_master_concat/'
src_name = 'PvP01_v1.filtered.core.coverage_filter_gtmissing.maf0.200.snps.for_FWS'

src_dir = '/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/vcf_master_concat/filter_core_coveragefilter_gtmissing_maf0.050/'
src_name = 'PvP01_v1.filtered.core.coverage_filter_gtmissing.maf0.050.snps'




vcf_file = paste0(src_dir, src_name, '.vcf.gz')
gds_file = paste0(src_dir, src_name, '.gds')

# # convert vcf to gds, and read in gds file
seqVCF2GDS(vcf_file, gds_file)
isolates <- seqOpen(gds_file)
fws <- getFws(isolates)


## some stats
hist(fws, breaks=100)
sum(fws < 0.95)
sum(fws >= 0.95)

## integrate with the masterfile with metadata
meta_data_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/data/PvGenomes_master.xlsx'
meta_data = read.xlsx(meta_data_file, sheet = 'PvGenomes_master_LatLong_filter')
meta_data = meta_data[-which(is.na(meta_data$Population)),]

meta_data$Population
head(meta_data)


# index_match = match(names(fws), meta_data$vcf_old)
# fws_meta_data = cbind.data.frame(meta_data[index_match,], fws)


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## make boxplot with fws data
ggplot(data=fws_meta_data, aes(x=Country, y=Fws)) + 
  geom_point(aes(color=Country), position="jitter") + 
  theme_bw() + 
  geom_hline(yintercept=0.95, linetype='dotted', col = 'black') +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1), 
        axis.text=element_text(size=14),
        legend.position = "none") 

ggplot(data=fws_meta_data, aes(x=Population, y=fws)) + 
  geom_point(aes(color=Population), position="jitter") + 
  theme_bw() + 
  geom_hline(yintercept=0.95, linetype='dotted', col = 'black') +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1), 
        axis.text=element_text(size=14),
        legend.position = "none") 



ggplot(data=meta_data, aes(x=Population, y=Fws, fill=Population)) + 
  # geom_boxplot(aes(color=Country)) + 
  geom_boxplot(color="black", lwd=.2) + 
  geom_point(data = meta_data[meta_data$Fws %in% boxplot.stats(meta_data$Fws)$out, ], 
             color = "black", size = 1, shape=21) +
  theme_bw() +
  geom_hline(yintercept=0.95, linetype='dotted', col = 'black') +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1), 
        axis.text=element_text(size=14),
        legend.position = "none") +
  scale_fill_manual(values = cbPalette) 

## Fws for only the South American countries
fws_meta_data_sam = fws_meta_data[fws_meta_data$Continent == "SAM",]
ggplot(data=fws_meta_data_sam, aes(x=Country, y=fws)) + 
  geom_point(aes(color=Country), position="jitter") + 
  theme_bw() + 
  geom_hline(yintercept=0.95, linetype='dotted', col = 'black') +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1), 
        axis.text=element_text(size=14),
        legend.position = "none") 

ggplot(data=fws_meta_data_sam, aes(x=Country, y=fws)) + 
  geom_boxplot(aes(color=Country)) + 
  theme_bw() +
  geom_hline(yintercept=0.95, linetype='dotted', col = 'black') +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1), 
        axis.text=element_text(size=14),
        legend.position = "none") 


## read in FWS data from the master file. This is after the Fws
## calculations have been run using moimix
meta_data_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/data/PvGenomes_master.xlsx'
# fws_meta_data = read.xlsx(meta_data_file, sheet = 'PvGenomes_master_LatLong_v2')
fws_meta_data = read.xlsx(meta_data_file, sheet = 'PvGenomes_master_LatLong_filter')
head(fws_meta_data)
fws_meta_data_sam = fws_meta_data[fws_meta_data$Continent == "SAM",]
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


ggplot(data=fws_meta_data_sam, aes(x=Country, y=Fws)) + 
  geom_point(aes(color=Country), position="jitter") + 
  theme_bw() + 
  geom_hline(yintercept=0.95, linetype='dotted', col = 'black') +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1), 
        axis.text=element_text(size=14),
        legend.position = "none") +
  scale_fill_manual(values = cbPalette) 

ggplot(data=fws_meta_data_sam, aes(x=Country, y=Fws, fill=Country)) + 
  # geom_boxplot(aes(color=Country)) + 
  geom_boxplot(color="black", lwd=.2) + 
  geom_point(data = fws_meta_data_sam[fws_meta_data_sam$Fws %in% boxplot.stats(fws_meta_data_sam$Fws)$out, ], 
             color = "black", size = 1, shape=21) +
  theme_bw() +
  geom_hline(yintercept=0.95, linetype='dotted', col = 'black') +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1), 
        axis.text=element_text(size=14),
        legend.position = "none") +
  scale_fill_manual(values = cbPalette) 
  
