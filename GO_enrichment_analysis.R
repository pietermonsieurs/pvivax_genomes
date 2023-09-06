library('openxlsx')

## set working directory
src_dir = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/results/GO/'
setwd(src_dir)

## define the different parameters to loop over
# databases = c('GO', 'KEGG')
databases = c('GO')


## define some parameters setting the minimum size of the GO category, as well
## as the minimum number of genes that should be detected in the genelist of 
## interest
min_go_size = 3
min_go_select_size = 2 

## specify the different input files
types = c('IBD_significance_E10', 'IBD_significance_E15')

for (database in databases) {
  
  ## create an empty list which will contain the different work sheets that
  ## can be written to an excel output file
  list_of_dfs = list()
  
  for (type in types) {
    
    ## create the GO file
    go_file = paste0(src_dir, type, '.GO.csv')
    
    ## read in the GO data sets
    go_data = read.table(go_file, sep=";", header=TRUE, quote="")
    head(go_data)
    dim(go_data)
    
    # filter on minimum size of the GO
    dim(go_data)
    sum(go_data$nr_in_GO >= min_go_size)
    go_data = go_data[go_data$nr_in_GO >= min_go_size,]
    go_data = go_data[go_data$nr_found >= min_go_select_size,]
    dim(go_data)
    
    go_data$pval = 1-phyper(go_data$nr_found, go_data$nr_sampled, go_data$genome_size-go_data$nr_sampled, go_data$nr_in_GO)
    go_data$pval = format(go_data$pval, scientific = FALSE)
    
    go_data[go_data$pval < 0.05,]
    dim(go_data[go_data$pval < 0.05,])
    
    ## export data to one .csv file per type of analysis. 
    csv_out = paste0(src_dir, type, '.GO.gsea.csv')
    
    write.table(go_data[go_data$pval < 0.05,], 
                sep=",", 
                file=csv_out)
    
    ## also write to a dictionarry so you can write everything in one
    ## go to an excel file
    sheet_name =  type
    list_of_dfs[[sheet_name]] = go_data[go_data$pval < 0.05,]   
    
  }
  
  ## write the output for a specific database (GO or KEGG) to an excel file
  ## containing different worksheets
  excel_out = paste0(database, ".enrichment.xlsx")
  write.xlsx(list_of_dfs, excel_out)
  
}
    



