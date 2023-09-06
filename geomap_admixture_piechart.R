library(ggmap)
library(maps)
library(countrycode)
library(rgeos)
library(rworldmap)
library(ggthemes)
library(ggrepel)
library(dplyr)
library(openxlsx)
library(scatterpie)


#### 1. create the world map #####

## file with the counts of samples per country. This has been manually 
## edited to correspond with the names in the official vocabulary of countries
country_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/data/2015_2020.Malaria_in_GeoSentinelDB.xlsx'
country_data = read.xlsx(country_file, sheet="combined")

# country_data = country_data[,-5] ## remove column with QC checks
country_data = country_data[-nrow(country_data),]
head(country_data)

## a complete list of all countries with their longitude and latitude
## which can be used to extract coordinates for the countries of interest
## has been downloaded from: https://developers.google.com/public-data/docs/canonical/countries_csv
coordinates_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/data/coordinates_countries_longitude_latitude.csv'
coord_data = read.csv(coordinates_file, sep=",")
coord_data[157,'country'] = "NA"
coord_data




## add the two letter abbreviation to get latitude and longitude in a second setup
## update of the countries need: China-Myanmar set to Myanmar and summed with the 
## other Myanmar samples (23 + 9 = 32), The philippines set to "Philippines" 
## and Papua Indonesia counted as Papua New Guinea
country_data$iso2c = countrycode(country_data$Country, "country.name", "iso2c")

## sort country data based on the abbreviation to match with the 
## longitude / altitude
country_data = country_data[order(country_data$iso2c),]

## DEBUG: check which country code not in coord_data
country_data[which(! country_data$iso2c %in% country_data$Country),]



## get long and lat coordinate based on the the iso2c code
country_data = cbind(country_data, coord_data[which(coord_data$country %in% country_data$iso2c),])


## map on the world 
world_data <- map_data("world")

# https://stackoverflow.com/questions/51350763/borders-and-colors-on-world-map-ggplot2
world <- map_data("world") 
world2 = merge(world, country_data, by.x="region", by.y="name", all.x=TRUE)


#ggplot(world2, aes(x = long, y = lat, fill=iso2c)) +
country_list = c(country_data$name, "Democratic Republic of the Congo", "Ivory Coast", "Myanmar", "Sao Tome and Principe", "South Sudan")
world_interest = world[world$region %in% country_list,]
country_data[which(! country_data$name %in% unique(world$region)),]

world_interest$interest=1


#### 2. create the country data - admixture group ####


# ---- 2.1. meta data ----
## read in the metadata, which is a copy of the pv genomes excel file
meta_data_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/data/PvGenomes_master.xlsx'
meta_data = read.xlsx(meta_data_file, sheet = "PvGenomes_master_LatLong_v2")
meta_data$loc_lang_lat = paste0(meta_data$Long, ":", meta_data$Lat)
head(meta_data)



# ---- 2.2 admixture data ----
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



# ---- 2.3 read in the admixture data ----
## read the assignment probabilities. This should be done
## using the K-value obtained from the cross validation 
## experiment 
K=2
K=3
K=4
K=5
K=8
K=10
K=11
admixture_dir = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/admixture/'
admixture_output_file = paste0(admixture_dir, 'PvP01_v1.filtered.core.coverage_filter_gtmissing.maf0.005.snps.ldpruned.filtered.',K,'.Q')
admixture_output_file = paste0(admixture_dir, 'PvP01_v1.filtered.core.coverage_filter_gtmissing.maf0.005.snps.ldpruned.filtered.',K,'.Q')
admixture_output_file = paste0(admixture_dir, 'PvP01_v1.filtered.core.coverage_filter_gtmissing.maf0.005.snps.SAM.ldpruned.filtered.',K,'.Q')

admixture_output = read.table(admixture_output_file)

## specify the column names, as they correspond to the probability that 
## a sample belongs to an admixture group
colnames_adm_output = paste0("K", seq(1,K,1))
colnames(admixture_output) = colnames_adm_output
rownames(admixture_output) = sample_order
# admixture_output = admixture_output[, -ncol(admixture_output)]
head(admixture_output)
dim(admixture_output)

# ---- 2.4 combinding the dataframe ----

## specify the dominant group as the cluster with the highes probability
## and add this as a separate column
admixture_output$group = apply(admixture_output,1,which.max) 
admixture_output$sample = sample_order
head(admixture_output)


## give rownames as sample names
## rownames(admixture_output) = meta_data_filter_sort$search_name_print # !!! WRONG !!!
index_for_meta_data = match(admixture_output$sample, meta_data_filter_sort$search_name_print)
plot_data = cbind(meta_data_filter_sort[index_for_meta_data, c('Population', 'Country', 'search_name_print', 'loc_lang_lat', 'year')], admixture_output)
sum(plot_data$search_name_print == plot_data$sample) == nrow(plot_data)
head(plot_data)
country_data = as.data.frame.matrix(table(plot_data$loc_lang_lat, plot_data$group))


# country_data = country_data[order(country_data$continent),]
# country_data$label = paste0(country_data$Country, "(", country_data$count, ")")
# country_data$color_order = seq(1, dim(country_data)[1])
# country_data$color_order = as.factor(country_data$color_order )
# country_data$GeoSentinel_increased = country_data$GeoSentinel + 20 # make dots more visual

## split the latitude and longitude from the 
coord_data = as.data.frame(str_split_fixed(rownames(country_data), ":", n=2))
colnames(coord_data) = c("longitude", "latitude")
coord_data$latitude = as.double(coord_data$latitude)
coord_data$longitude = as.double(coord_data$longitude)

country_data_plot = cbind(country_data, coord_data)

## columns to take into account depend on the number of Ks
columns_to_include = seq(1,K,1)
country_data_plot$count=rowSums(country_data_plot[,columns_to_include])
country_data_plot = country_data_plot[-nrow(country_data_plot),]


cols = brewer.pal(n = 8, name = "Dark2")
# cols = palette.colors(n = 9, palette = "Okabe-Ito") # https://jfly.uni-koeln.de/color/#pallet
col.palette <- colorRampPalette(cols,
                                space = "Lab")

ggplot() +   
  theme(panel.background = element_rect('light blue')) + 
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill="gray70", color="gray50", size=0.4, alpha=0.20) +
  geom_polygon(data = world_interest, aes(x = long, y = lat, group = group), color = "gray50", fill="gray20", size=0.4, alpha=0.20) +
  theme_map() +
  coord_cartesian(ylim=c(-30,20), xlim=c(-110,-40)) + 
  geom_scatterpie(data = country_data_plot, aes(x = longitude, y = latitude, r=log10(count+5/3)),
                  cols=columns_to_include,
                  color="gray30", alpha=0.90) + 
  scale_fill_manual(values = col.palette(K+1)) 
  


# + 
#   scale_fill_manual(values = col.palette(K+1))

  
  # geom_scatterpie(data = country_data, aes(x = longitude, y = latitude, r=log10(count+5/3)), 
  # geom_scatterpie(data = country_data, aes(x = longitude, y = latitude), 
                  # cols=c('GeoSentinel', 'Malariagen'),
                  cols=columns_to_include,
                  color="gray30", alpha=0.90) #  +
# geom_scatterpie_legend(country_data$count, x=-100, y=-20) #+

# geom_label_repel(data=arrange(country_data, country) ,aes(x=longitude, y=latitude, label=label, color=color_order), 
#                 show.legend=FALSE, force=50,force_pull = 10,
#                 min.segment.length=0,
#                 #fontface="bold",
#                 label.size=NA,
#                 fill=NA,
#                 alpha=1,
#                 seed=1234,
#                 label.padding=0.1,
#                 max.overlaps = Inf
#                 
#               ) +
#geom_point( data=country_data, aes(x=longitude, y=latitude, size=Malariagen, color=color_order)) +
# geom_point( data=country_data, aes(x=longitude, y=latitude, size=Malariagen), color="black", pch=21) +

# geom_point( data=country_data, aes(x=longitude, y=latitude, size=non_fp, color=color_order)) +
# geom_point( data=country_data, aes(x=longitude, y=latitude, size=non_fp), color="black", pch=21) + 
# 
# theme(legend.position="none")



country_data






## only one colour code - no colouring per country

ggplot() +   
  theme(panel.background = element_rect('light blue')) + 
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill="gray95", color="gray60", size=0.4, alpha=0.20) +
  geom_polygon(data = world_interest, aes(x = long, y = lat, group = group), color = "gray70", fill="gray20", size=0.4, alpha=0.20) +
  theme_map() +
  coord_cartesian(ylim=c(-30,40), xlim=c(-110,140)) + 
  #geom_text(data=country_data,aes(x=longitude, y=latitude, size=count, label=name) ) + 
  # 
  # geom_label_repel(data=arrange(country_data, country) ,aes(x=longitude, y=latitude, label=label), 
  #                  show.legend=FALSE, force=50,force_pull = 10,
  #                  min.segment.length=0,
  #                  #fontface="bold",
  #                  label.size=NA,
  #                  alpha=0.70,
  #                  seed=1234,
  #                  label.padding=0
  # )  +

# geom_label_repel(data=arrange(country_data, country) ,aes(x=longitude, y=latitude, label=label),
#                 show.legend=FALSE, force=50,force_pull = 10,
#                   min.segment.length=0,
#                   #fontface="bold",
#                   label.size=NA,
#                   fill=NA,
#                   alpha=1,
#                   seed=1234,
#                   label.padding=0,
#                   max.overlaps = Inf
# ) +
geom_point( data=country_data, aes(x=longitude, y=latitude, size=GeoSentinel), color="dodgerblue4", alpha=0.70)



## check with pieplots: https://guangchuangyu.github.io/2016/12/scatterpie-for-plotting-pies-on-ggplot/




## add color code based on the amount of samples per country

## add to world_interest the number of samples
world_interest$count = 0
world_interest$count = country_data$GeoSentinel[match(world_interest$region, country_data$name)]
world_interest$count_malariagen = country_data$Malariagen[match(world_interest$region, country_data$name)]



## add to world_malariagen the countries with data in the 
## malariagen Pf6 database
world_interest$malariagen = NA
world_interest[world_interest$region %in% country_data$Country[country_data$Malariagen > 0],]$malariagen = 1



ggplot() +   
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill="gray90", color="gray60", size=0.4, alpha=0.80) +
  geom_polygon(data = world_interest, aes(x = long, y = lat, group = group, fill=log10(count_malariagen)), alpha=0.80, color = "gray70", size=0.4) +
  #geom_polygon(data = world_interest, aes(x = long, y = lat, group = group, fill=malariagen), alpha=0.80, color = "gray80", size=0.4) +
  
  # scale_fill_distiller(palette = "YlOrRd", trans = "reverse") +
  scale_fill_distiller(palette = "Blues", trans = "reverse") +
  
  theme_map() +
  coord_cartesian(ylim=c(-30,40), xlim=c(-110,140)) +
  
  # geom_point( data=country_data, aes(x=longitude, y=latitude, size=log10(GeoSentinel+ 10), color=color_order)) +
  # geom_point( data=country_data, aes(x=longitude, y=latitude, size=log10(GeoSentinel+ 10)), color="brown", alpha=.5) +
  # geom_point( data=country_data, aes(x=longitude, y=latitude, size=log10(GeoSentinel + 10)), color="black", pch=21) #  +
  
  geom_point( data=country_data, aes(x=longitude, y=latitude, size=GeoSentinel_increased), color="brown", alpha=.5) +
  geom_point( data=country_data, aes(x=longitude, y=latitude, size=GeoSentinel_increased), color="black", pch=21)+ 
  
  guides(size=guide_legend("GeoSentinel"), fill = "none")
# guides(size=guide_legend("GeoSentinel"), fill = guide_legend("MalariaGen"))


# geom_point( data=country_data, aes(x=longitude, y=latitude, size=non_fp, color=color_order)) +
# geom_point( data=country_data, aes(x=longitude, y=latitude, size=non_fp), color="black", pch=21) + 
# 


theme(legend.position="none")
#geom_text(data=country_data,aes(x=longitude, y=latitude, size=count, label=name) ) + 

# geom_label_repel(data=country_data, aes(x=longitude, y=latitude, label=label), 
#                  show.legend=FALSE,
#                  min.segment.length=0.50,
#                  # fontface="bold",
#                  label.size=NA,
#                  alpha=0.70,
#                  seed=1234,
#                  label.padding=0,
#                  point.padding=0) +

# geom_label_repel(data=country_data ,aes(x=longitude, y=latitude, label=label), 
#                  show.legend=FALSE,
#                  min.segment.length=0.50,
#                  # fontface="bold",
#                  label.size=NA,
#                  fill=NA,
#                  alpha=1, 
#                  seed=1234,
#                  label.padding=0,
#                  point.padding=0,
#                  max.overlaps=Inf)


#geom_point( data=country_data, aes(x=longitude, y=latitude, size=count, alpha=0.70)) 



#### make pie charts 



