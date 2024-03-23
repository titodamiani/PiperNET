# Usage: python blastp.py -db <path_to_blast_database> -q <path_to_query_sequence> --mcov 70 --msim 50 --proteome <path_to_proteome> -o <output_fasta>
#   -db, --blastDB: Path to the BLAST database.
#   -q, --query: Path to the query sequence (FASTA file).
#   --proteome: Path to the proteome sequence database (FASTA file).
#   --mcov: Minimum coverage for a BLAST hit to be retained (default: 0).
#   --msim: Minimum similarity for a BLAST hit to be retained (default: 0).
#   -o, --output: Output folder where matching sequences will be stored in a FASTA file.


import argparse
import pandas as pd
from Bio import SeqIO
from utils.rnaseq import blastp

def main():
    parser = argparse.ArgumentParser(description="XXX")
    parser.add_argument("-db", "--blastDB", dest="blastDB", required=True, help="Path to BLAST database.")    
    parser.add_argument("-q", "--query", dest="query", required=True, help="Path to query serquence (FASTA file).")
    parser.add_argument("--proteome", dest="proteome", required=True, help="Path to proteome (FASTA file).")
    parser.add_argument("--mcov", dest="min_cov", default=0, help="Minimum coverage for a BLAST hit to be retained.")
    parser.add_argument("--msim", dest="min_sim", default=0, help="Minimum similarity for a BLAST hit to be retained.")
    parser.add_argument("-o", "--output", dest="output_fasta", required=True, help="Output folder")
    args = parser.parse_args()

    #run blastp
    blast_results = blastp(args.blastDB, args.query, min_cov = float(args.min_cov), min_sim = float(args.min_sim))
    print(blast_results)
    hit_ids = set(blast_results['seq_id'])

    #empty list to store the matching sequences
    hit_sequences = []

    #retreive hit sequences from proteome.pep
    with open(args.proteome, 'r') as proteome_handle:
        for record in SeqIO.parse(proteome_handle, 'fasta'):
            if record.id in hit_ids:
                hit_sequences.append(record)

    # Write the matching sequences to the output FASTA file
    with open(args.output_fasta, 'w') as output_handle:
        for record in hit_sequences:
            full_id = blast_results.loc[blast_results['seq_id'] == record.id, 'full_id'].values[0]
            print(full_id)
            record.id = full_id
            record.description = ""  #clear record description field
            SeqIO.write(record, output_handle, 'fasta-2line')
            

    print(f"{len(hit_sequences)} BLAST hits stored to '{args.output_fasta}'.")


if __name__ == "__main__":
    main()