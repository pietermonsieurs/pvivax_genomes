#!/usr/bin/env python3

from Bio import Entrez

def get_biosample_accession(run_accession):
    try:
        Entrez.email = "pmonsieurs@itg.be"  # Replace with your email
        print(f"run_accession =={run_accession}==")
        handle = Entrez.esummary(db="sra", id=run_accession)
        record = Entrez.read(handle)
        handle.close()
        return record[0]['Runs'][0]['sra_biosample']
    except Exception as e:
        print(f"Error retrieving BioSample for {run_accession}: {e}")
        return None

def main():
    # run_accessions = [  # Replace with your list of run accession numbers
    #     "ERR5740735",
    #     "SRR12134582",
    #     "SRR12134581",
    #     "SRR12134579"] 
    run_accessions = [  # Replace with your list of run accession numbers
        "SRR000001", "SRR000002", "SRR000003",  # Add more accession numbers as needed
    ]

    biosample_accessions = {}
    for run_accession in run_accessions:
        biosample_accession = get_biosample_accession(run_accession)
        if biosample_accession:
            biosample_accessions[run_accession] = biosample_accession

    for run_accession, biosample_accession in biosample_accessions.items():
        print(f"Run Accession: {run_accession}, BioSample Accession: {biosample_accession}")

if __name__ == "__main__":
    main()
