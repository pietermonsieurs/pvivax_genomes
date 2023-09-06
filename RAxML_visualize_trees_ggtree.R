library(ggtree)
library(treeio)
library(tidyr)
library(dplyr)
library(naniar)
library(openxlsx)
library(tidyverse)

# if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install()
# 
# ## biocLite("BiocUpgrade") ## you may need this
# BiocManager::install("ggtree")


## input parameters
data_dir = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/raxml/'
setwd(data_dir)

## tree for the population structure worldwide
tree <- read.newick('PvP01_v1.filtered.maf_0.200.snps.core.plink_filtered.T3.raxml.bestTree')

## tree for the population structure SAM
tree <- read.newick('PvP01_v1.filtered.core.coverage_filter_gtmissing.maf0.005.snps.SAM.T3.raxml.bestTree')


tree

## create simple ggtree
p = ggtree(tree)
p = p + theme_tree2()
p = p + geom_tiplab(size=2)
p


## convert to tree to a tibble data frame which can then 
## be used to add additional information
x <- as_tibble(tree)


## read in the master file with the meta data
country_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/data/PvGenomes_master.xlsx'
country_data = read.xlsx(country_file, sheet='PvGenomes_master_LatLong_v2')
country_data

## for some of the Vietnam samples, the sample name is very complex, including
## a slash. Only keep the part after the slash
# head(pca_data)
# problematic_indices = grep("/", pca_data[,1])
# problematic_names = pca_data[problematic_indices,1]
# 
# correct_names <- sapply(strsplit(problematic_names, "/"), function(x) x[length(x)])
# pca_data[problematic_indices, 1] = correct_names

x$country = country_data[match(x$label, country_data$search_name_print),]$Country
x$continent = country_data[match(x$label, country_data$search_name_print),]$Continent
x$population = country_data[match(x$label, country_data$search_name_print),]$Population
x


## join the metadata with the intial tree
x_new = full_join(tree,x)



## convert the tibble tree back to a phylogenetic tree
## that can be visualized using ggtree

# tree = as.phylo(x_new)


## coloring for worldwide scale
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pops = sort(unique(x$population))
pop_colors = cbPalette[1:length(pops)]
pop_colors = setNames(pop_colors, pops)

p = ggtree(x_new)
p = p + geom_tree(linewidth = .01) 
p = p + geom_tippoint(aes(color=population), size=.5)
# p = p + geom_tippoint(aes(color=country), size=.5)
p = p + scale_color_manual(values = pop_colors)  
p = p + guides(color = guide_legend(override.aes = list(size = 5)))
p = p + theme(legend.text=element_text(size=14))
p

### save the output for the world tree
output_dir = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/pictures_pub/'
output_file = paste0(output_dir, 'tree_world.png')
output_file_rds = paste0(output_dir, 'tree_world.rds')
ggsave(output_file, plot=p, width=12, height = 3, dpi=300)
saveRDS(p, output_file_rds)





## coloring for SAM scale
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pops = sort(unique(x$country))
pop_colors = cbPalette[1:length(pops)]
pop_colors = setNames(pop_colors, pops)

p = ggtree(x_new)
p = p + geom_tree(linewidth = .01) 
# p = p + geom_tippoint(aes(color=population), size=.5)
p = p + geom_tippoint(aes(color=country), size=.5)
p = p + scale_color_manual(values = pop_colors)  
p = p + guides(color = guide_legend(override.aes = list(size = 5)))
p = p + theme(legend.text=element_text(size=14))
p

### save the output for the world tree
output_dir = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/pictures_pub/'
output_file = paste0(output_dir, 'tree_sam.png')
output_file_rds = paste0(output_dir, 'tree_sam.rds')
ggsave(output_file, plot=p, width=6, height = 3, dpi=300)
saveRDS(p, output_file_rds)


#### some tryout snippets form previous experiments

# p = ggtree(x_new, layout='circular')
# p = layout(p, layout = "circular", circular.radius = 0.8)
# p = p + geom_tiplab(aes(angle = angle, filter = isTip), size = 0.01, hjust =-100) 
p  
# p = p + theme_tree2()
# p = p + geom_tiplab(size=2)
# p = p + geom_tippoint(mapping = aes(color=continent))
# p = p + geom_tippoint(mapping = aes(color=population))
# p = p + geom_tippoint(mapping = aes(color=country))
# p = p + geom_tippoint(aes(color=population), size=.5)
# p = p + scale_color_manual(values = pop_colors)  
#   
# p = p + guides(color = guide_legend(override.aes = list(size = 5)))
# p = p + theme(legend.text=element_text(size=8))
# p = p + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
#               axis.line = element_blank(),
#         axis.text.x = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks = element_blank())
# p = p + geom_tree(linewidth = 0.5) 
# 
# # p = p + xlim(-10, NA)
# 
# p
# 
# p <- layout(p, layout = "circular", node.width = 0, node.height = 0, tip.width = 0, tip.height = 0)


## focus on the West Asia countries. Keep the population for
## all samples, except for the west-asia samples where you print
## the country, to make the outgroup of the WAS samples visible 
## in the tree
# x = as.data.frame(x)
x$focus = NA
x[which(x$population == "WAS"),]$focus = x[which(x$population == "WAS"),]$country
x_new = full_join(tree,x)


# p = ggtree(x_new, layout='circular')
p = ggtree(x_new)
p = p + theme_tree2()
# p = p + geom_tiplab(size=2)
# p = p + geom_tippoint(mapping = aes(color=continent))
# p = p + geom_tippoint(mapping = aes(color=population))
p = p + geom_tippoint(mapping = aes(color=focus), size=0.5)
p = p + guides(color = guide_legend(override.aes = list(size = 14)))
p = p + theme(legend.text=element_text(size=14))
# p = p + geom_tiplab2(size=2.5, color="black", align=TRUE)

p


