#!/usr/bin/env python3

## gff file containing the size of the chromosomes
gff_file = '/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/refgenomes/pvivax/PlasmoDB-46_PvivaxP01.gff'
region_list = '/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/gvcf_master_merged/chrom_list_regions.csv'

## create the chromosome list
chrom_list = list()
for chrom in range(1,15):
    chrom_list.append(f"{chrom:02d}")

chrom_list.append('MIT')
chrom_list.append('API')
print(chrom_list)

## open the output file and print the header containing the different
## column headers
out_fh = open(region_list, 'w')
header = "chrom,start,end\n"
out_fh.write(header)

## define the window size in nucleotides: this defines how big the 
## different chunks of data are that will be used as input
## for the genotyping
window_size = 100000

## loop over the different chromosomes and look up the size
## of the chromosome in the .gff file
for chrom in chrom_list:
    
    ## define the full chromosome name
    full_chrom = f'PvP01_{chrom}_v1'

    ## find the size of the chromosome
    gff_fh = open(gff_file, 'r')
    for line in gff_fh:
        line = line.rstrip()
        if line.startswith(f"##sequence-region {full_chrom} "):
            chrom_size = int(line.split(" ")[3])
            print(f"chrom {full_chrom} -> {chrom_size}") 

    ## split into different windows
    end_of_chrom_reached = 0
    start = 1
    while not end_of_chrom_reached:
        end = start + window_size - 1
        if end > chrom_size:
            end = chrom_size
            end_of_chrom_reached = 1
        
        out_line = f"{chrom},{start},{end}\n"
        out_fh.write(out_line)

        start = start + window_size 

out_fh.close()
