export PATH=/user/antwerpen/205/vsc20587/data/software/python_lib/bin:${PATH}
export PYTHONPATH=/user/antwerpen/205/vsc20587/data/software/python_lib/lib/python3.8/site-packages:${PYTHONPATH}
module load Python/3 


## go to data fastq directory and get the accession numbers
cd /user/antwerpen/205/vsc20587/scratch/plasmodium_pvgenomes/data/fastq_valdivia/
while read -r accnr
do 
    echo ${accnr}
    ffq $accnr -o ${accnr}.json
done < sra_numbers.csv
