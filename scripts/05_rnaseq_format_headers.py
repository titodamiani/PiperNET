# Usage: python format_headers.py -i <input_proteome> -o <output_proteome> --prefix <header_prefix>
#   -i, --input: input rnaSPAdes proteome from transXpress (transcriptome.pep)
#   --prefix: Prefix to add to each header after '>'.

from pathlib import Path
import pandas as pd
from src.rnaseq_utils import get_config
import argparse

def main():
    parser = argparse.ArgumentParser(description="Format rnaSPAdes haeders.")
    parser.add_argument("--config", dest='config_path', help="Path to config.yaml", default='config/rnaseq_dataprep.yaml')
    args = parser.parse_args()


    #load sample folder paths from config.yaml
    assembly_paths = get_config(args.config_path, data='rna-seq', file='assembly') #dictionary
    prefixes = get_config(args.config_path, data='rna-seq', file='prefix')


    #iterate over path dictionary
    for input, prefix in zip(assembly_paths.values(), prefixes.values()):

        try:
            
            #i/o file paths
            input = input.with_name('transcriptome_clstr.pep')
            output = input.with_name('proteome.pep')

            print(f'Formatting headers in {input}...')
            with open(input, 'r') as handle, open(output, 'w') as output_csv:
                
                #initialize header and sequence
                header = None
                sequence = ''
                
                for line in handle:

                    #header line
                    if line.startswith('>'):

                        #if not the first header, write previous sequence to output file
                        if header is not None:
                            output_csv.write(sequence + '\n')
                            sequence = ''  # Reset the sequence

                        # Store the header for later use
                        header = line

                        # Get the first part of the header before the first white space
                        header = header.split()[0]
                        header = header.split('_')[6:]
                        header = '_'.join(header[0:])

                        # Add the specified header prefix
                        new_header = '>' + str(prefix) + header

                        #write new header to output file
                        output_csv.write(new_header + '\n')

                    else:
                        #remove leading/trailing spaces from sequence line
                        line = line.strip()

                        #concatenate sequence lines into a single string
                        sequence += line

                # Write the last sequence
                if header is not None:
                    output_csv.write(sequence + '\n')
                
            print(f'Formatting complete! Formatted proteome written to {output}')

        #skip to next file if input is not found
        except FileNotFoundError:
            print(f'{input} not found. Skipping to next sample.')
            continue

if __name__ == "__main__":
    main()