## split the full url file into batches of 100 url

cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/fastq_malariagen
source_file=download_url.csv

split -l 100 -d --additional-suffix=.txt $source_file batch

## remove header url from the first batch
vi batch00.txt

## run per batch 
while read url; do wget $url; done < batch00.txt
while read url; do wget $url; done < batch01.txt
while read url; do wget $url; done < batch02.txt
while read url; do wget $url; done < batch03.txt
while read url; do wget $url; done < batch04.txt
while read url; do wget $url; done < batch05.txt
while read url; do wget $url; done < batch06.txt
while read url; do wget $url; done < batch07.txt
while read url; do wget $url; done < batch08.txt
while read url; do wget $url; done < batch09.txt
while read url; do wget $url; done < batch10.txt
while read url; do wget $url; done < batch11.txt

## run additionally one batch with the missing accession numbers
## that were not retrieved in extraction of URL using SRA-explorerd
## This list can be generated using the script 
## inputdata_download_find_missing_malariagen.py
while read url; do wget $url; done < download_url_missed.csv

## manually create ftp donwload URL for ERR2299660
while read url; do wget $url; done < download_url_missed_ERR2299660.csv