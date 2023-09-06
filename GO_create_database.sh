
## download the .gaf database from PlasmoDB
cd /Users/pmonsieurs/programming/data/GO
curl -L  https://plasmodb.org/common/downloads/Current_Release/PvivaxP01/gaf/PlasmoDB-64_PvivaxP01_GO.gaf -O PlasmoDB-64_PvivaxP01_GO.gaf
curl -L  https://plasmodb.org/common/downloads/Current_Release/PvivaxP01/gaf/PlasmoDB-64_PvivaxP01_Curated_GO.gaf -O PlasmoDB-64_PvivaxP01_Curated_GO.gaf

## convert the .gaf to .gmt version using a script already
## used for leishmania. Data are stored in data/GO/. run for 
## both normal and curated gaf file
/Users/pmonsieurs/programming/plasmodium_pvgenomes/bin/GO_parse_gaf.py 


## run the script to go from a file containg a list of genes
## to a text file that can be used to calculate the hypergeometric
## distribution
cd /Users/pmonsieurs/programming/plasmodium_pvgenomes/bin

./GO_genes_to_GOenrichment.py -g /Users/pmonsieurs/programming/plasmodium_pvgenomes/results/GO/IBD_significance_E10.csv -d ~/programming/data/GO/PlasmoDB-64_PvivaxP01_GO.gmt -o /Users/pmonsieurs/programming/plasmodium_pvgenomes/results/GO/IBD_significance_E10.GO.csv

./GO_genes_to_GOenrichment.py -g /Users/pmonsieurs/programming/plasmodium_pvgenomes/results/GO/IBD_significance_E15.csv -d ~/programming/data/GO/PlasmoDB-64_PvivaxP01_GO.gmt -o /Users/pmonsieurs/programming/plasmodium_pvgenomes/results/GO/IBD_significance_E15.GO.csv