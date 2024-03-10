import argparse
from src.rnaseq_utils import *
from pathlib import Path
import pandas as pd
import re
import csv
from Bio import SeqIO
from io import StringIO


def main():
    parser = argparse.ArgumentParser(description="Parse transcriptome.pep and merge transXpress annotations.")
    parser.add_argument("--config", dest='config_path', help="Path to config.yaml", default='config/rnaseq_dataprep.yaml')
    args = parser.parse_args()

    #load sample folder paths from config.yaml
    assembly_paths = get_config(args.config_path, data='assembly')
    prefixes = get_config(args.config_path, data='prefix', filepaths=False)
    
    #iterate over path dictionary
    for input, prefix in zip(assembly_paths.values(), prefixes.values()):
        
        try:
            #i/o file paths        
            input = input.with_name('transcriptome_clstr.pep')
            output = input.with_name('proteome.csv')
            pfam_path = input.parent / 'annotations' / 'pfam_orfs.out'
            signalp_path = input.parent / 'annotations' / 'signalp_pep.out'
            blastp_path = input.parent / 'annotations' /'sprotblastp_orfs.out'
            targetp_path = input.parent / 'annotations' / 'targetp_pep.out'
            tmhmm_path = input.parent / 'annotations' / 'tmhmm_pep.out'
            kallisto_path = input.parent / 'transcriptome_expression_isoform.tsv'


            #parse transcriptome_clstr.pep
            print(f'Parsing {input}...')
            proteome = []

            with open(input) as handle:
                for record in SeqIO.parse(handle, "fasta"):

                    header = str(record.description) #full rnaSPAdes header
                    sequence = str(record.seq)
                    id = header.split('_', 6)[6].split(' ', 1)[0] #split after 6th '_' and before whitespace
                    
                    #length and coverage of original transcript
                    header_splits = header.split('_')
                    transcript_length = header_splits[3]
                    transcript_coverage = header_splits[5]

                    #orf info
                    orf_type = extract_orf_type(header)
                    orientation = extract_orientation(header)
                    orf_score = extract_score(header)
                    
                    new_row = {
                        'id': id,
                        'sequence': sequence,
                        'length': len(sequence), #protein sequence length
                        'orf_type': orf_type,
                        'orientation': orientation,
                        'orf_score': orf_score,
                        'transcript_length': transcript_length,
                        'transcript_coverage': transcript_coverage
                    }
                    proteome.append(new_row)

                proteome = pd.DataFrame(proteome)

            #parse Pfam
            pfam={}

            with open(pfam_path) as handle:
                for line in handle:
                    
                    if (line.startswith("#")): continue
                    row = re.split(" +", line, 22)
                    if row[3] not in pfam:
                        pfam[row[3]] = [row[1], row[22].rstrip('\n'), row[6]]
                

            pfam = pd.DataFrame.from_dict(pfam, orient='index', columns=['pfam', 'domain', 'pfam_e-value']).reset_index().rename(columns={'index': 'id'}) #convert to dataframe
            pfam = shorten_id(pfam)



            #parse BLASTp
            blastp={}

            with open(blastp_path) as handle:
                csv_reader = csv.reader(handle, delimiter="\t")
                for row in csv_reader:
                    if (len(row) < 13): continue
                    blastp[row[0]] = [row[12], row[10]]

            blastp = pd.DataFrame.from_dict(blastp, orient='index', columns=['blastp', 'blastp_e-value']).reset_index().rename(columns={'index': 'id'}) #convert to dataframe
            blastp = shorten_id(blastp)

            #parse SignalP
            with open(signalp_path, 'r') as handle:
                signalp = ''.join(line for line in handle if not line.startswith('#'))

            signalp = pd.read_csv(StringIO(signalp), sep='\s+', usecols=[0, 1], names=["id", "signalp"])
            signalp = shorten_id(signalp)

            #parse TargetP
            with open(targetp_path, 'r') as handle:
                targetp_data = ''.join(line for line in handle if not line.startswith('#'))

            targetp = pd.read_csv(StringIO(targetp_data), sep='\s+', usecols=[0, 1], names=["id", "targetp"])
            targetp = shorten_id(targetp)

            #parse TMHMM
            tmhmm = pd.read_csv(tmhmm_path, sep='\t', header=None)
            tmhmm.columns = ['id', 'tmhmm', 'topology'] #rename the columns
            tmhmm['tmhmm'] = tmhmm['tmhmm'].str.replace('PredHel=', '') #remove 'PredHel='
            tmhmm['topology'] = tmhmm['topology'].str.replace('Topology=', '') #remove 'Topology='
            tmhmm = shorten_id(tmhmm)


            #merge annotations
            dfs = [pfam, blastp, signalp, targetp, tmhmm]
            for df in dfs:
                proteome = proteome.merge(df, on='id', how='left')
            
            #merge tpm (kallisto)
            kallisto = pd.read_csv(kallisto_path, sep='\t')
            kallisto.rename(columns={'Unnamed: 0':'id'}, inplace=True)
            kallisto = shorten_id(kallisto).set_index('id') #set id as index
            proteome['temp_id'] = proteome['id'].apply(lambda x: x.split('.')[0]) #create temporary_id w/o ORF
            proteome = proteome.merge(kallisto, how='inner', left_on='temp_id', right_index=True) #merge expression
            proteome.drop(['temp_id'], axis=1, inplace=True)
            

            #add prefix to id column
            proteome['id'] = prefix + proteome['id']
            
            #write to csv
            proteome.to_csv(output, index=False)
            print(f'Parsing complete! Proteome written to {output}')


        #skip to next file if input is not found
        except FileNotFoundError:
            print(f'{input} not found. Skipping to next sample.')
            continue

if __name__ == "__main__":
    main()