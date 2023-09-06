#!/usr/bin/env python3

import pandas as pd

master_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/data/PvGenomes_master.xlsx'
master_file_out = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/data/PvGenomes_master.with_geo.xlsx'
malariagen_file = '/Users/pmonsieurs/programming/plasmodium_pvgenomes/data/metadata_2022_malariagen.xlsx'

## read in the master data file
master_data = pd.read_excel(master_file)
print(master_data.head)

## read in the malariagen file
malaria_data = pd.read_excel(malariagen_file)
print(malaria_data.head)

## do simple merge 
master_data_with_geo = master_data.merge(malaria_data[['Sample', 'Lat', 'Long']], 
    how='left', 
    left_on='identifier', 
    right_on='Sample')
print(master_data_with_geo.head)
master_data_with_geo.to_excel(master_file_out)





