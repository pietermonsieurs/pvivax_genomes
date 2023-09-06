library(ggplot2)
library(scatterpie)
library(dplyr)
library(RColorBrewer)


#### 1. create input data set: metadata + admixture ####


# ---- 1.1. meta data ----
## read in the metadata, which is a copy of the pv genomes excel file
meta_data_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/data/PvGenomes_master.xlsx'
meta_data = read.xlsx(meta_data_file, sheet = "PvGenomes_master_LatLong_filter")
meta_data$loc_lang_lat = paste0(meta_data$Long, ":", meta_data$Lat)
head(meta_data)

# ---- 1.2 admixture data ----
## the .fam file specifies the order of the samples in the output of the 
## admixture output. The last one is the admixture file create for the 
## South America plot
sample_order_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/admixture/PvP01_v1.filtered.maf_0.005.snps.core.ldpruned.filtered.fam'
sample_order_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/admixture/PvP01_v1.filtered.core.coverage_filter_gtmissing.maf0.005.snps.ldpruned.filtered.fam'
sample_order_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/admixture/PvP01_v1.filtered.core.coverage_filter_gtmissing.maf0.005.snps.SAM.ldpruned.filtered.fam'
sample_order = read.table(sample_order_file, sep="\t")[,2]
sample_order


## filter + order the metadata
meta_data_filter = meta_data[meta_data$search_name_print %in% sample_order,]
meta_data_filter_sort= meta_data_filter[order(meta_data_filter$search_name_print),]
dim(meta_data_filter_sort)


head(meta_data_filter_sort)



# ---- 1.3 read in the admixture data ----
## read the assignment probabilities. This should be done
## using the K-value obtained from the cross validation 
## experiment 
K=11
admixture_dir = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/admixture/'
# admixture_output_file = paste0(admixture_dir, 'PvP01_v1.filtered.core.coverage_filter_gtmissing.maf0.005.snps.ldpruned.filtered.',K,'.Q')
# admixture_output_file = paste0(admixture_dir, 'PvP01_v1.filtered.core.coverage_filter_gtmissing.maf0.005.snps.ldpruned.filtered.',K,'.Q')
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

# ---- 1.4 combinding the dataframe ----

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


#### 2. create Piechart bubble plot ####

# ---- 2.1 conversion to input format ----

plot_data_src = plot_data
plot_data = plot_data_src


## explore the data to see which years should be combined
table(plot_data$year)
plot_data[is.na(plot_data$year),]
plot_data[plot_data$year < 2007,]$year = "< 2007"
plot_data$year = as.factor(plot_data$year)

## create a new data frame that can be used to store all 
## the data for the piechart plot
# pie_data = cbind.data.frame(unique(plot_data$year), unique(plot_data$Country))
pie_data = expand.grid(unique(plot_data$year), unique(plot_data$Country))
colnames(pie_data) = c('year', 'country')

pie_data

## add additional columsn that convert the year and the 
## country into a numeric value. However, this sorts the year
## in alphabetical order...
pie_data <- pie_data %>% 
 mutate(year_num = as.numeric(as.factor(year)), 
        country_num = as.numeric(as.factor(country)))

pie_data <- pie_data %>% 
  mutate(country_num = as.numeric(as.factor(country)))


## option 2: add data frame that contains a conversion 
## table between year and numeric value
year_conversion_df = cbind.data.frame(levels(as.factor(plot_data$year)), seq(1,length(levels(as.factor(plot_data$year)))))
colnames(year_conversion_df) = c('year', 'year_num')
pie_data$year_num = year_conversion_df[match(pie_data$year, year_conversion_df$year),]$year_num


## create a new data format
# pie_data = plot_data[, c('Country', 'year', paste0("K", seq(1,11,1)), 'group')]
pie_data <- cbind(pie_data, K1 = 0, K2 = 0, K3 = 0, K4 = 0, K5 = 0, K6 = 0, K7 = 0, K8 = 0, K9 = 0, K10 = 0, K11 = 0)
for (i in 1:nrow(plot_data)) {
  year = plot_data[i,'year']
  country = plot_data[i,'Country']
  admix_group = plot_data[i,'group']
  admix_column = paste0("K", admix_group)
  print(paste0("before: ", pie_data[pie_data$year == year & pie_data$country == country, admix_column] ))
  pie_data[pie_data$year == year & pie_data$country == country, admix_column] = pie_data[pie_data$year == year & pie_data$country == country, admix_column] + 1
  print(paste0("after: ", pie_data[pie_data$year == year & pie_data$country == country, admix_column] ))
}
pie_data


# ---- 2.2 piechart bubble plot ----

##  https://stackoverflow.com/questions/66535303/geom-scatterpie-with-non-numeric-axes
columns_for_pie = paste0("K", seq(1,11,1))
cols_breaks = seq(1,length(unique(pie_data$year)))
cols_labels = sort(levels(as.factor(pie_data$year)))
rows_breaks = seq(1,length(unique(pie_data$country)))
rows_labels = unique(pie_data$country)

## specify colors 
cols = brewer.pal(n = 8, name = "Dark2")
col.palette <- colorRampPalette(cols,
                                space = "Lab")

## add the rowSums to get the radius of the pie charts
pie_data$size = rowSums(pie_data[,5:15])

p <- ggplot() + 
  geom_scatterpie(data = pie_data, aes(x=year_num, y=country_num, r=log(size+1, base=1000)),
                  cols=columns_for_pie, alpha=0.80) + 
  geom_scatterpie_legend(log(pie_data$size+2, base=1000), x=15, y=1, labeller = function(x) (as.integer(1000 ** (x)))) +
  scale_fill_manual(values=col.palette(12)) + 
  scale_x_continuous(labels = cols_labels, breaks = cols_breaks) + 
  scale_y_continuous(labels = rows_labels, breaks = rows_breaks) + 
  xlab('year') + ylab('country') + 
  theme_bw() + 
  theme(
    axis.text = element_text(size = 14),        # Increase axis text size
    axis.title = element_text(size = 14),       # Increase axis title size
    legend.text = element_text(size = 14),      # Increase legend text size
    legend.title = element_text(size = 14))



p


pie_data$country


set.seed(123)
long <- rnorm(50, sd=100)
lat <- rnorm(50, sd=50)
d <- data.frame(long=long, lat=lat)
d <- with(d, d[abs(long) < 150 & abs(lat) < 70,])
n <- nrow(d)
d$region <- factor(1:n)
d$A <- abs(rnorm(n, sd=1))
d$B <- abs(rnorm(n, sd=2))
d$C <- abs(rnorm(n, sd=3))
d$D <- abs(rnorm(n, sd=4))
d[1, 4:7] <- d[1, 4:7] * 3
head(d)

d$region[1:20] = 1
d$lat[1:20] = d$lat[1]
d$long[1:20] = d$long[1]


ggplot() + geom_scatterpie(aes(x=long, y=lat, group=region), data=d,
                           cols=LETTERS[1:4]) + coord_equal()



