library(ggplot2)
library(cowplot)
library(ggpubr)
library(dplyr)
library(reshape)
library(openxlsx)
library(RColorBrewer)

## install package specifically for plotting admixture analysis, which 
## can be loaded using Netview
# install.packages("remotes")
# remotes::install_github("esteinig/netview")
# library(netview)
# 
# install.packages('popkin')
# library(popkin)
# 
# plot_admix(admixture_output_file,
#            col = RColorBrewer::brewer.pal(max(ncol(Q), 3), "Paired"))
# 
# http://pophelper.com/ -> maximum size of 200kb => size of Q20 is too high (250kb)


#### 1. input data ####

# ---- 1.1. meta data ----
## read in the metadata, which is a copy of the pv genomes excel file
meta_data_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/data/PvGenomes_master.xlsx'
meta_data = read.xlsx(meta_data_file, sheet = "PvGenomes_master_LatLong_v2")
head(meta_data)


# ---- 1.2 admixture data ----
## the .fam file specifies the order of the samples in the output of the 
## admixture output 
sample_order_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/admixture/PvP01_v1.filtered.maf_0.005.snps.core.ldpruned.filtered.fam'
sample_order_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/admixture/PvP01_v1.filtered.core.coverage_filter_gtmissing.maf0.005.snps.ldpruned.filtered.fam'
sample_order_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/admixture/PvP01_v1.filtered.core.coverage_filter_gtmissing.maf0.005.snps.SAM.ldpruned.filtered.fam'
sample_order = read.table(sample_order_file, sep="\t")[,2]
sample_order


## filter + order the metadata
meta_data_filter = meta_data[meta_data$search_name_print %in% sample_order,]
meta_data_filter_sort= meta_data_filter[order(meta_data_filter$search_name_print),]


# sum(sample_order== meta_data_filter_sort$search_name_print)
# length(sample_order)
# sample_order[919:923]
# meta_data_filter_sort[919:923,]$search_name_print
# 
# 
# sample_order[919:930]
# meta_data_filter_sort[919:930,]$search_name_print


# ---- 1.3 read in the admixture data ----
## read the assignment probabilities. This should be done
## using the K-value obtained from the cross validation 
## experiment 
K=2
K=3
K=5
K=11
K=10
K=21
admixture_dir = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/admixture/'
admixture_output_file = paste0(admixture_dir, 'PvP01_v1.filtered.core.coverage_filter_gtmissing.maf0.005.snps.ldpruned.filtered.', K, '.Q')
admixture_output_file = paste0(admixture_dir, 'PvP01_v1.filtered.core.coverage_filter_gtmissing.maf0.005.snps.SAM.ldpruned.filtered.', K, '.Q')

admixture_output = read.table(admixture_output_file)

## specify the column names, as they correspond to the probability that 
## a sample belongs to an admixture group. Also add the rownames by 
## the sample_order
colnames_adm_output = paste0("K", seq(1,K,1))
colnames(admixture_output) = colnames_adm_output
rownames(admixture_output) = sample_order
admixture_output$sample = sample_order
# admixture_output = admixture_output[, -ncol(admixture_output)]
head(admixture_output)
dim(admixture_output)


## using NetView
# plotAdmixture(admixture_output_file, graph=NULL, structurePlot=FALSE, palette="Dark2", colourN=8)
# plotAdmixture(admixture_output_file, metaData=meta_data_filter)


# d <- dist(admixture_output)
# h <- hclust(d)
# dend <- as.dendrogram(h)
# plot(dend)


#### 2. visualisation ####

# ---- 2.1 sorting the samples according to the groups ----

## specify the dominant group as the cluster with the highes probability
## and add this as a separate column
admixture_output$group = apply(admixture_output,1,which.max) 

## give rownames as sample names
##rownames(admixture_output) = meta_data_filter_sort$search_name_print # !!! WRONG !!!
index_for_meta_data = match(admixture_output$sample, meta_data_filter_sort$search_name_print)
plot_data = cbind(meta_data_filter_sort[index_for_meta_data, c('Population', 'Country', 'search_name_print', 'year')], admixture_output)
sum(plot_data$search_name_print == plot_data$sample) == nrow(plot_data)
head(plot_data)




## sort the different samples according to the group where they belong to 
## and next to the probability per group: ugly solution - still find a way 
## to sort by multiple columsn names which are defined as variables
# if (K==5) {plot_data_sorted = plot_data[with(plot_data, order(group, -K1, -K2, -K3, -K4, -K5)),]}
# if (K==10) {plot_data_sorted = plot_data[with(plot_data, order(group, -K1, -K2, -K3, -K4, -K5, -K6, -K7, -K8, -K9, -K10)),]}
# if (K==11) {plot_data_sorted = plot_data[with(plot_data, order(group, -K1, -K2, -K3, -K4, -K5, -K6, -K7, -K8, -K9, -K10, -K11)),]}
# if (K==20) {plot_data_sorted = plot_data[with(plot_data, order(group, -K1, -K2, -K3, -K4, -K5, -K6, -K7, -K8, -K9, -K10, -K11, -K12, -K13, -K14, -K15, -K16, -K17, -K18, -K19, -K20)),]}
# if (K==21) {plot_data_sorted = plot_data[with(plot_data, order(group, -K1, -K2, -K3, -K4, -K5, -K6, -K7, -K8, -K9, -K10, -K11, -K12, -K13, -K14, -K15, -K16, -K17, -K18, -K19, -K20, -K21)),]}


## alternative approach for plot_data_sorted, where you select per group
## the subset of data, and sort them according to the group that is selected. 
## this gives nicer plots than the code above, as you sort per dominant group
## according to that group
plot_data_sorted = data.frame()
for (K_value in 1:K) {
  plot_data_sub = plot_data[plot_data$group == K_value,]
  colname = paste0("K", K_value)
  plot_data_sub_sorted = plot_data_sub[order(plot_data_sub[,colname], decreasing = TRUE),]
  if (K_value == 1) {
    plot_data_sorted = plot_data_sub_sorted
  }else{
    plot_data_sorted = rbind.data.frame(plot_data_sorted, plot_data_sub_sorted)
  }
}
head(plot_data_sorted)


## remove the last column (group) as otherwise this will be
## added to the plot data
plot_data_sorted = plot_data_sorted[,-ncol(plot_data_sorted)]
head(plot_data_sorted)


# ---- 2.2 plot the different groups in a barplot ----

## melt the plot data to work with ggplot
plot_data_sorted_melt = melt(plot_data_sorted, id=c('Population', 'Country', 'search_name_print', 'year', 'sample'))
plot_data_sorted_melt$search_name_print = factor(plot_data_sorted_melt$search_name_print, levels=plot_data_sorted$search_name_print)
colnames(plot_data_sorted_melt)[dim(plot_data_sorted_melt)[2]] = 'Group'
head(plot_data_sorted_melt)

## specify the colors. One way is to work with a set of predefined
## colors and use e.g. scale_fill_brewer(palette="Set3"), but not always
## enough colors present. So you can create your own color scheme starting from 
## exsiting one and use "scale_fill_manual(values = ...)
cols = brewer.pal(n = 8, name = "Dark2")
# cols = palette.colors(n = 9, palette = "Okabe-Ito") # https://jfly.uni-koeln.de/color/#pallet
col.palette <- colorRampPalette(cols,
                                space = "Lab")


## create the admixture plot
p_admix = ggplot(plot_data_sorted_melt, aes(fill=variable, y=Group, x=search_name_print)) + 
    geom_bar(position="fill", stat="identity") + 
    theme(axis.text.x=element_blank(), #remove x axis labels
          axis.ticks.x=element_blank()) + 
    # scale_fill_manual(values = col.palette(K+1))
    scale_fill_manual(values = col.palette(K+1)) + 
    xlab('Sample') + ylab('Ancestry')
  

    # scale_fill_brewer(palette="Set3")
p_admix


# ---- 2.3 create the side bars ----

## create the side bar for the bar plot to show the country or the 
plot_data_artificial = plot_data_sorted
plot_data_artificial$search_name_print = factor(plot_data_artificial$search_name_print, levels=plot_data_artificial$search_name_print)


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot_data_artificial$dummy = 1
p_population = ggplot(plot_data_artificial, aes(fill=Population, y=dummy, x=search_name_print)) + 
  geom_bar(position="fill", stat="identity") + 
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + 
  scale_fill_manual(values = cbPalette) +
  xlab(NULL) + ylab(NULL)
p_population
  #scale_fill_brewer(palette="Set3")
  

p_country = ggplot(plot_data_artificial, aes(fill=Country, y=dummy, x=search_name_print)) + 
  geom_bar(position="fill", stat="identity") + 
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + 
        scale_fill_manual(values = cbPalette) + 
  xlab(NULL) + ylab(NULL)
p_country

p_year = ggplot(plot_data_artificial, aes(fill=year, y=dummy, x=search_name_print)) + 
  geom_bar(position="fill", stat="identity") + 
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank()) + 
  scale_fill_continuous(limits=c(2000,2020), 
                        breaks=seq(2000,2020,10),
                        
                        guide=guide_colourbar(reverse = TRUE))
  # scale_fill_gradientn(limits = c(2000, 2020),
  #                       breaks=seq(2000,2020,1),
  #                      colours=c("grey90", viridis::inferno(10000, begin = 0, end = 1, direction = -1)))
                        # colours = rev(heat.colors(length(seq(2000,2020,1)))))


# ---- 2.4 combine different plots ----
p_side = p_country
# p_side = p_population

ggarrange(p_side, p_admix, ncol = 1, nrow = 2, heights = c(1, 10))

# Combine the two plots using the plot_grid() function from cowplot
combined_plots <- plot_grid(p_side + theme(legend.position = "none"), p_admix + theme(legend.position = "none"), ncol = 1, align = "v", rel_heights = c(0.2, 0.8))
combined_plots

# Create a legend for the combined plot
side_legend <- get_legend(p_side)
admix_legend <- get_legend(p_admix)

# Add the legend to the combined plot and place it on the right side
final_plot_admixture <- plot_grid(combined_plots, admix_legend, side_legend, ncol = 3, rel_widths = c(0.82, 0.08, 0.10), align = "h", axis='tb')

# Display the final plot
final_plot_admixture
output_dir = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/pictures_pub/'

## specify the output for SAM
output_file = paste0(output_dir, 'admixture_SAM.png')
output_file_rds = paste0(output_dir, 'admixture_SAM.rds')
ggsave(output_file, plot=final_plot_admixture, width=12, height = 3, dpi=300)
saveRDS(final_plot_admixture, output_file_rds)

## specify the output for world
output_file = paste0(output_dir, 'admixture_world.png')
output_file_rds = paste0(output_dir, 'admixture_world.rds')
ggsave(output_file, plot=final_plot_admixture, width=12, height = 3, dpi=300)
saveRDS(final_plot_admixture, output_file_rds)




#### 3. PCA visualisation ####

# ---- 3.1 read in the data and link metadata ----
# pca_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/admixture/pca_results.SAM.eigenvec'
pca_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/admixture/PvP01_v1.no_outliers.pca_results.eigenvec' # SAM
pca_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/admixture/pca_results.eigenvec' # World

pca_data = read.table(pca_file, sep="\t")

## for the world database, there is an additinoal column in the 
## plink output file, and problematic naming of one of the samples from
## Mauritania
if (pca_file == '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/admixture/pca_results.eigenvec') {
  pca_data = pca_data[,-1]
  pca_data = pca_data[- which(pca_data[,1] == "Mauritania1"),]
}

## for some of the Vietnam samples, the sample name is very complex, including
## a slash. Only keep the part after the slash
head(pca_data)
problematic_indices = grep("/", pca_data[,1])
problematic_names = pca_data[problematic_indices,1]

correct_names <- sapply(strsplit(problematic_names, "/"), function(x) x[length(x)])
pca_data[problematic_indices, 1] = correct_names


# head(meta_data_filter_sort)
# head(cbind(pca_data, meta_data_filter_sort))

## select only those sample which are present in the meta data and
## sort them according to the sample name (search_name_print)
meta_data_filter = meta_data[meta_data$search_name_print %in% pca_data[,1],] ## for SAM
# meta_data_filter = meta_data[meta_data$search_name_print %in% pca_data[,1] ## for world
meta_data

meta_data_filter_sort= meta_data_filter[order(meta_data_filter$search_name_print),]
head(meta_data_filter_sort)

pca_data[! pca_data[,1] %in% meta_data$search_name_print,]


## also sort the PCA data according to the sample name and sort 
## them according to first column (=sample name). 
pca_data_sorted = pca_data[order(pca_data[,1]),]

## Also check whether the sample names are in the same order
## before merging
sum(pca_data_sorted[,1] == meta_data_filter_sort$search_name_print)
dim(pca_data_sorted)

## do combination of both data frames
plot_data_pca = cbind(pca_data_sorted, meta_data_filter_sort)
head(plot_data_pca)
colnames(plot_data_pca)[2:6] = c('PC1', 'PC2', 'PC3', 'PC4', 'PC5')
colnames(plot_data_pca)[2:6] = c('PC1', 'PC2', 'PC3', 'PC4', 'PC5')
head(plot_data_pca)
dim(plot_data_pca)



# ---- 3.2 make plot using Country level analysis (SAM) ---

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pca_plot_PC1_2_SAM = ggplot(plot_data_pca, aes(x=PC1,y=PC2)) + 
    geom_point(aes(fill=Country), alpha=1, colour="black",pch=21, size=2) + 
    #geom_text(aes(label=V1), vjust=1,hjust=1) +
    scale_fill_manual(values = cbPalette) + 
    theme_bw()
pca_plot_PC1_2_SAM

ggplot(plot_data_pca, aes(x=PC1,y=PC2)) + 
  geom_point(aes(fill=year), alpha=1, colour="black",pch=21, size=3) + 
  #geom_text(aes(label=V1), vjust=1,hjust=1) +
  theme_bw() + 
  scale_fill_continuous(limits=c(2000,2020), 
                        breaks=seq(2000,2020,10),
                        guide=guide_colourbar(reverse = TRUE))

ggplot(plot_data_pca, aes(x=PC3,y=PC4)) + 
  geom_point(aes(fill=Country), alpha=1, colour="black",pch=21, size=3) + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw()


## specify the output for SAM
output_dir = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/pictures_pub/'
output_file = paste0(output_dir, 'pca_1_2_SAM.png')
output_file_rds = paste0(output_dir, 'pca_1_2_SAM.rds')
ggsave(output_file, plot=pca_plot_PC1_2_SAM, width=12, height = 3, dpi=300)
saveRDS(pca_plot_PC1_2_SAM, output_file_rds)



# ---- 3.3 make plot using Population level analysis ---

pca_plot_PC1_2_world = ggplot(plot_data_pca, aes(x=PC1,y=PC2)) + 
  geom_point(aes(fill=Population), alpha=1, colour="black",pch=21, size=2) + 
  #geom_text(aes(label=V1), vjust=1,hjust=1) +
  scale_fill_manual(values = cbPalette) + 
  theme_bw()
pca_plot_PC1_2_world

ggplot(plot_data_pca, aes(x=PC1,y=PC2)) + 
  geom_point(aes(fill=year), alpha=1, colour="black",pch=21, size=3) + 
  #geom_text(aes(label=V1), vjust=1,hjust=1) +
  theme_bw() + 
  scale_fill_continuous(limits=c(2000,2020), 
                        breaks=seq(2000,2020,10),
                        guide=guide_colourbar(reverse = TRUE))

ggplot(plot_data_pca, aes(x=PC3,y=PC4)) + 
  geom_point(aes(fill=Country), alpha=1, colour="black",pch=21, size=3) + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw()


## specify the output for world
output_dir = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/pictures_pub/'
output_file = paste0(output_dir, 'pca_1_2_world.png')
output_file_rds = paste0(output_dir, 'pca_1_2_world.rds')
ggsave(output_file, plot=pca_plot_PC1_2_world, width=12, height = 3, dpi=300)
saveRDS(pca_plot_PC1_2_world, output_file_rds)




# --- 3.4 add admixture information to PCA plot? ----

## to be implemented ##



## functions to be loaded? ###
plot.admixture<-function(directory){
  
  require(reshape2)
  require(grid)
  require(ggplot2)
  require(plyr)
  require(RColorBrewer)
  
  Qmatrix <- dir(directory,pattern="*.Q")
  fam<-dir(directory, pattern="*.fam")
  temp.name<-read.csv(paste(directory, fam, sep=""), sep="\t", header=FALSE)
  
  imax<-function(x,n){
    tmp<-as.data.frame(cbind(x, 1:length(x)))
    colnames(tmp)<-c("val","index")
    stmp<-tmp[order(tmp$val, decreasing=TRUE),]
    return(stmp$index[n])
  }
  vmax<-function(x,n){
    tmpd<-sort(x, decreasing=TRUE)
    return(tmpd[n])
  }
  
  parse.Q.files<-function(x){
    dat<-read.csv(file=paste(directory,x, sep=""), sep=" ",header=FALSE)
    ldat<-length(dat)
    
    outter <- NULL	
    
    for(i in 1:ldat){	
      outter<-c(outter,list(apply(dat, 1, FUN=imax, n=i)))
      outter<-c(outter,list(apply(dat, 1, FUN=vmax, n=i)))
    }
    
    dat<-cbind(temp.name$V2,dat)
    colnames(dat)<-c( "names", 1:ldat)
    
    dat<-dat[do.call(order, outter),]
    dat$ov<-1:length(dat$names)
    
    dat2<-melt(dat, id.vars=c("names", "ov"))	
    
    
    dat2$Krun<-ldat
    colnames(dat2)<-c("Name","ov","Admixture.Group","Value", "Krun")
    
    
    return(dat2)
  }	
  
  
  
  plotgg<-function(datframe){
    
    my.max.k<-max(as.numeric(as.character(datframe$Admixture.Group)))
    
    
    datframe<-datframe[order(datframe$ov),]
    datframe$Name<-factor(datframe$Name, levels=datframe$Name, ordered=TRUE)
    
    my.col<-colors()[c(26,547,498,69,33,51,536,100,76,200,300,400,450)]
    
    the.plot<-ggplot(datframe, aes(x=Name, y=Value, fill=Admixture.Group))+geom_bar(stat="identity")
    the.plot<-the.plot+scale_fill_manual(values = my.col[1:my.max.k], name="Admixture group")
    the.plot<-the.plot+theme_classic(18)
    the.plot<-the.plot+theme(axis.text.x = element_text(angle = 90, hjust = 1))
    the.plot<-the.plot+labs(x="Individual", y="fraction ancestry")
    the.plot
    
  }
  
  fdat<-ldply(Qmatrix, parse.Q.files, .inform=TRUE)
  
  my.plots<-dlply(fdat, .(Krun), plotgg)
  
  return(my.plots)
  
}



plot.all.together<-function(file.name, plotCols, list.of.plots){
  
  numPlots = length(list.of.plots)
  plotRows = ceiling(numPlots/plotCols)
  # Fiddle with the to adjust your plot dimentions 
  pdf(file=paste(file.name, "pdf", sep="."),bg="transparent", width=18*plotCols, height=8*plotRows)
  
  grid.newpage() 
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols))) 
  vplayout <- function(x, y) 
    viewport(layout.pos.row = x, layout.pos.col = y) 
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(my.plots[[i]], vp = vplayout(curRow, curCol ))
  }
  
  dev.off() 
  
}

