#!/usr/bin/env python3

import sys

src_dir = '/Users/pmonsieurs/programming/data/GO/'

## Ldon
# gaf_file_in = src_dir + 'TriTrypDB-48_LdonovaniBPK282A1_GO.gaf'
# gtm_file_out = src_dir + 'TriTrypDB-48_LdonovaniBPK282A1_GO.gtm'

## Lbraz
# gaf_file_in = src_dir + 'TriTrypDB-55_LbraziliensisMHOMBR75M2904_2019_GO.gaf'
# gtm_file_out = src_dir + 'TriTrypDB-55_LbraziliensisMHOMBR75M2904_2019_GO.gtm'

## Pv curated
# gaf_file_in = f"{src_dir}/PlasmoDB-64_PvivaxP01_Curated_GO.gaf"
# gtm_file_out = f"{src_dir}/PlasmoDB-64_PvivaxP01_Curated_GO.gmt"

## Pv all
gaf_file_in = f"{src_dir}/PlasmoDB-64_PvivaxP01_GO.gaf"
gtm_file_out = f"{src_dir}/PlasmoDB-64_PvivaxP01_GO.gmt"

go_description_file = src_dir + 'go.obo'

go_full = open(go_description_file, 'r')
go_names = {}
for line in go_full:
    if line.startswith("id: GO:"):
        line = line.rstrip()
        go_class = line.replace("id: ", "")
        name_line = next(go_full)
        name_line = name_line.rstrip()
        go_name = name_line.replace("name: ", "")
        go_names[go_class] = go_name
        # print(go_class)

# print(go_names)
# sys.exit()


gaf_in = open(gaf_file_in, 'r')
gtm_out = open(gtm_file_out, 'w')
go = {}
for line in gaf_in:
    if line.startswith("!"): continue
    data = line.split("\t")
    GO_class = data[4]
    gene = data[1].replace(".1", "")
    print(GO_class)

    if GO_class in go.keys():
        go[GO_class].append(gene)
    else:
        go[GO_class] = [gene]

# print(go.keys())

for go_class in go.keys():
    genes = "\t".join(go[go_class])
    if go_class not in go_names.keys():
        print(f"go_class {go_class} not found in GO names")
        # sys.exit()
        continue
    go_name = go_names[go_class]
    line_out = f"{go_class}\t{go_name}\t{genes}\n"
    print(line_out)
    gtm_out.write(line_out)


