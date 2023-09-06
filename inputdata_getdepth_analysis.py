#!/usr/bin/env python3

import os
import sys
import pandas as pd


### option 1: read all depth files from the depth directory ###
## however, this takes a lot of time, and might be easier / faster
## to read from. Might be faster to integrate into the slurm command

# depth_dir = '/user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/results/depth_from_gvcf/'

# for depth_file in os.listdir(depth_dir):
#     depth_file = f"{depth_dir}/{depth_file}"
#     print(depth_file)
#     depth_data = pd.read_csv(depth_file)
#     # depth_data.iloc[:,2] = depth_data.iloc[:,2].replace(".", 0)
#     depth_data.iloc[:,2] = pd.to_numeric(depth_data.iloc[:,2], errors='coerce')
#     print(depth_data.iloc[:,2].median())
#     print(sum(depth_data.iloc[:,2]>=5))
#     print(depth_data.shape)


### option 2: read depth file from command line ###

## minimum coverage to be "callable"
min_cov = 5

## read input and define output
depth_file = sys.argv[1]
depth_file_stats = depth_file.replace(".csv", ".stats.csv")
print(depth_file_stats)

## open output_file
out_fh = open(depth_file_stats, 'w')


## read in the data and convert to integers so that 
## values indicated as "." are converted to NaN
depth_data = pd.read_csv(depth_file)
# depth_data.iloc[:,2] = depth_data.iloc[:,2].replace(".", 0)
depth_data.iloc[:,2] = pd.to_numeric(depth_data.iloc[:,2], errors='coerce')

## get median and percentage callable
median_depth = depth_data.iloc[:,2].median()
pos_above_cutoff = sum(depth_data.iloc[:,2] >= min_cov)
pos_all = depth_data.shape[0]
perc_callable = 100*pos_above_cutoff/pos_all

## write out
sample = os.path.basename(depth_file)
line_out = f"{sample},{median_depth},{pos_above_cutoff},{pos_all},{perc_callable}\n"
print(line_out)
out_fh.write(line_out)

## close filehandle
out_fh.close()
