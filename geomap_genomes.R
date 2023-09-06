library(ggmap)
library(maps)
library(countrycode)
library(rgeos)
library(rworldmap)
library(ggthemes)
library(ggrepel)
library(dplyr)
library(openxlsx)



## file with the counts of samples per country. This has been manually 
## edited to correspond with the names in the official vocabulary of countries
# country_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/data/pv_genomes_input_ggmap.csv'
# country_data = read.csv(country_file, sep=",")
country_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/data/PvGenomes_master.xlsx'
country_data = read.xlsx(country_file, sheet='summary_per_country_geomap')
country_data = country_data[-nrow(country_data),] ## remove "total" field
country_data



## a complete list of all countries with their longitude and latitude
## which can be used to extract coordinates for the countries of interest
## has been downloaded from: https://developers.google.com/public-data/docs/canonical/countries_csv
coordinates_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/data/coordinates_countries_longitude_latitude.csv'
coord_data = read.csv(coordinates_file, sep=",")
coord_data


## add the two letter abbreviation to get latitude and longitude in a second setup
## update of the countries need: China-Myanmar set to Myanmar and summed with the 
## other Myanmar samples (23 + 9 = 32), The philippines set to "Philippines" 
## and Papua Indonesia counted as Papua New Guinea. 
# country_data$iso2c = countrycode(country_data$Country, "country.name", "iso2c")
country_data$iso2c = countrycode(country_data$country, "country.name", "iso2c")
country_data = country_data[! is.na(country_data$iso2c),]

## sort country data based on the abbreviation to match with the 
## longitude / altitude
country_data = country_data[order(country_data$iso2c),]

## get long and lat coordinate based on the the iso2c code
country_data = cbind(country_data, coord_data[which(coord_data$country %in% country_data$iso2c),])

## previous step integrates a new "country" column, but this is a copy 
## of the iso2c country code, so remove from the data to prevent error
## when plotting

country_data = country_data[,-7]


## map on the world 
world_data <- map_data("world")
# p = ggplot() + geom_polygon(data = world_data, aes(x=long, y = lat, group = group), fill="grey", alpha=0.8)
# p = 
# p = p + theme_void() 
# p
# 
# p = ggplot() + geom_polygon(data = world_data, aes(x=long, y = lat, group = group),
#                             fill="grey", alpha=0.5, border="black")
# p = p + theme_void() 
# p = p + coord_cartesian(ylim=c(-30,30))
# p = p + geom_point( data=country_data, aes(x=longitude, y=latitude, size=count)) 
#   
# p
# 
# 
# mapWorld <- borders("world", colour="black", fill="gray80", size=0.2) 
# mp <- ggplot() +   mapWorld
# mp = mp + geom_point( data=country_data, aes(x=longitude, y=latitude, size=count)) 
# mp = mp + coord_cartesian(ylim=c(-30,30))
# mp = mp + theme_bw()
# mp


# https://stackoverflow.com/questions/51350763/borders-and-colors-on-world-map-ggplot2
world <- map_data("world") 
world2 = merge(world, country_data, by.x="region", by.y="name", all.x=TRUE)

# 
# ggplot(world2, aes(x = long, y = lat, fill=iso2c)) +
#   geom_polygon(col = "white") + 
#   theme_map() +
#   coord_cartesian(ylim=c(-30,50))
#   

#ggplot(world2, aes(x = long, y = lat, fill=iso2c)) +
world_interest = world[world$region %in% country_data$name,]
world_interest$interest=1
  

country_data = country_data[order(country_data$continent),]
country_data$label = paste0(country_data$Country, "(", country_data$count, ")")
country_data$color_order = seq(1, dim(country_data)[1])
country_data$color_order = as.factor(country_data$color_order )

## update for countries with in-house data
for (i in 1:nrow(country_data)) {
  print(country_data$inhouse[i])
  if (country_data$inhouse[i] > 0) {
    print("updating label")
    country_data$label[i] = paste0(country_data$Country[i], "(", country_data$count[i], ")")
    # country_data$label[i] = paste0(country_data$Country[i], "(", 
    #                                country_data$count[i] - country_data$inhouse[i], 
    #                                "+",country_data$inhouse[i], ")")
    # #print(paste0(country_data$Country[i], "(", 
                 #country_data[i]$count - country_data[i]$inhouse, 
                 #" + ",country_data[i]$inhouse, ")")
    #)
  }
  
}





ggplot() +   
  theme(panel.background = element_rect('light blue')) + 
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill="gray95", color="gray60", size=0.4, alpha=0.20) +
  geom_polygon(data = world_interest, aes(x = long, y = lat, group = group), color = "gray70", fill="gray20", size=0.4, alpha=0.20) +
  theme_map() +
  coord_cartesian(ylim=c(-30,40)) + 
  #geom_text(data=country_data,aes(x=longitude, y=latitude, size=count, label=name) ) + 
  
  geom_label_repel(data=arrange(country_data, continent) ,aes(x=longitude, y=latitude, label=label, color=color_order), 
  # geom_label_repel(data=country_data,aes(x=longitude, y=latitude, label=label, color=color_order), 
                   show.legend=FALSE, force=50,force_pull = 10,
                   min.segment.length=0,
                   #fontface="bold",
                   label.size=NA,
                   alpha=0.70,
                   seed=1234,
                   label.padding=0.1
  ) +  
  
  geom_label_repel(data=arrange(country_data, continent) ,aes(x=longitude, y=latitude, label=label, color=color_order), 
                  show.legend=FALSE, force=50,force_pull = 10,
                  min.segment.length=0,
                  #fontface="bold",
                  label.size=NA,
                  fill=NA,
                  alpha=1,
                  seed=1234,
                  label.padding=0.1
                  
                ) +
  geom_point( data=country_data, aes(x=longitude, y=latitude, size=count, color=color_order)) +
  geom_point( data=country_data, aes(x=longitude, y=latitude, size=count), color="black", pch=21) + 
  theme(legend.position="none")

  

country_data






## only one colour code - no colouring per country

ggplot() +   
  theme(panel.background = element_rect('light blue')) + 
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill="gray95", color="gray60", size=0.4, alpha=0.20) +
  geom_polygon(data = world_interest, aes(x = long, y = lat, group = group), color = "gray70", fill="gray20", size=0.4, alpha=0.20) +
  theme_map() +
  coord_cartesian(ylim=c(-30,40)) + 
  #geom_text(data=country_data,aes(x=longitude, y=latitude, size=count, label=name) ) + 
  
  geom_label_repel(data=arrange(country_data, continent) ,aes(x=longitude, y=latitude, label=label), 
                   show.legend=FALSE, force=50,force_pull = 10,
                   min.segment.length=0,
                   #fontface="bold",
                   label.size=NA,
                   alpha=0.70,
                   seed=1234,
                   label.padding=0
  )  +
  
  geom_label_repel(data=arrange(country_data, continent) ,aes(x=longitude, y=latitude, label=label), 
                   show.legend=FALSE, force=50,force_pull = 10,
                   min.segment.length=0,
                   #fontface="bold",
                   label.size=NA,
                   fill=NA,
                   alpha=1,
                   seed=1234,
                   label.padding=0
                   
  ) +
  geom_point( data=country_data, aes(x=longitude, y=latitude, size=count, alpha=0.70)) 




## add color code based on the amount of samples per country

## add to world_interest the number of samples
world_interest$count = 0
world_interest$count = country_data$count[match(world_interest$region, country_data$name)]


ggplot() +   
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill="gray95", color="gray60", size=0.4, alpha=0.20) +
  geom_polygon(data = world_interest, aes(x = long, y = lat, group = group, fill=log10(count)), alpha=0.80, color = "gray70", size=0.4) +
  # scale_fill_distiller(palette = "YlOrRd", trans = "reverse") +
  scale_fill_distiller(palette = "Blues", trans = "reverse") +
  
  theme_map() +
  coord_cartesian(ylim=c(-30,40)) + 
  #geom_text(data=country_data,aes(x=longitude, y=latitude, size=count, label=name) ) + 
  
  geom_label_repel(data=country_data, aes(x=longitude, y=latitude, label=label), 
                   show.legend=FALSE,
                   min.segment.length=0.50,
                   # fontface="bold",
                   label.size=NA,
                   alpha=0.70,
                   seed=1234,
                   label.padding=0,
                   point.padding=0) +
  
  geom_label_repel(data=country_data ,aes(x=longitude, y=latitude, label=label), 
                   show.legend=FALSE,
                   min.segment.length=0.50,
                   # fontface="bold",
                   label.size=NA,
                   fill=NA,
                   alpha=1, 
                   seed=1234,
                   label.padding=0,
                   point.padding=0)

                   
  
  #geom_point( data=country_data, aes(x=longitude, y=latitude, size=count, alpha=0.70)) 




#### fine grained analysis of South America #####

## read in the raw data from excel
country_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/data/PvGenomes_master.xlsx'
country_data = read.xlsx(country_file, sheet='PvGenomes_master_LatLong')
# country_data = country_data[-nrow(country_data),] ## remove "total" field
country_data

## concatenate the latitude and longitude for identification of the location
country_data$geo = paste0(country_data$Lat, ":", country_data$Long)
country_data_sam = country_data[country_data$Continent == 'SAM',]
country_data_sam

plot_data = as.data.frame(table(country_data_sam$geo))
colnames(plot_data) = c('geo', 'count')
plot_data

plot_data$Country = NA
plot_data$Lat = NA
plot_data$Long = NA
for (i in 1:nrow(plot_data)) {
  geo = plot_data[i, 'geo']
  meta_data = country_data_sam[country_data_sam$geo == geo,][1,]
  plot_data[i, 'Lat'] = as.double(meta_data$Lat[1])
  plot_data[i, 'Long'] = as.double(meta_data$Long[1])
  plot_data[i, 'Country'] = meta_data$Country[1]
}

plot_data



ggplot() +   
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill="gray95", color="gray60", size=0.4, alpha=0.20) +
  # geom_polygon(data = world_interest, aes(x = long, y = lat, group = group, fill=log10(count)), alpha=0.80, color = "gray70", size=0.4) +
  # scale_fill_distiller(palette = "YlOrRd", trans = "reverse") +
  scale_fill_distiller(palette = "Blues", trans = "reverse") +
  
  theme_map() +
  coord_cartesian(ylim=c(-20,20), xlim=c(-100,-40)) + 
  #geom_text(data=country_data,aes(x=longitude, y=latitude, size=count, label=name) ) + 
  
  geom_label_repel(data=plot_data, aes(x=Long, y=Lat, label=count), 
                   show.legend=FALSE,
                   min.segment.length=0.50,
                   # fontface="bold",
                   label.size=NA,
                   alpha=0.70,
                   seed=1234,
                   label.padding=0,
                   point.padding=0) +
  
  geom_label_repel(data=plot_data ,aes(x=Long, y=Lat, label=count), 
                   show.legend=FALSE,
                   min.segment.length=0.50,
                   # fontface="bold",
                   label.size=NA,
                   fill=NA,
                   alpha=1, 
                   seed=1234,
                   label.padding=0,
                   point.padding=0) + 
  geom_point( data=plot_data, aes(x=Long, y=Lat, size=count, color=Country), alpha=0.70)



