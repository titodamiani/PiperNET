import subprocess
import argparse
from src.rnaseq_utils import get_config
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Create BLAST databases.")
    parser.add_argument("--config", dest='config_path', help="Path to config.yaml", default='config/rnaseq_dataprep.yaml')
    args = parser.parse_args()

    #load sample folder paths from config.yaml
    assembly_paths = get_config(args.config_path, data='proteome_fasta') #dictionary
    species_names = get_config(args.config_path, data='species', filepaths=False)

    #iterate over path dictionary
    for input, species in zip(assembly_paths.values(), species_names.values()):
        
        try:
            #i/o file paths
            # input = input.with_name('proteome.pep')
            blastDB_dir = input.parent / 'blastDB'
            output = blastDB_dir / input.parent.stem

            # check if input file exists
            if not input.exists():
                raise FileNotFoundError()

            #create blastDB directory if not exist
            if not blastDB_dir.exists():
                blastDB_dir.mkdir(parents=True, exist_ok=True)
                print(f"Directory '{blastDB_dir}' created.")

            else:
                print(f"Directory '{blastDB_dir}' already exists.")


            #construct command
            command = f'makeblastdb -in {input} -dbtype prot -out {output} -title {species} -parse_seqids'
            print(f'BLAST database created for {species} at {blastDB_dir}')
            
            #run command
            subprocess.run(command, shell=True)

        #skip to next file if input is not found
        except FileNotFoundError:
            print(f'{input} not found. Skipping to next sample.')
            continue

if __name__ == "__main__":
    main()