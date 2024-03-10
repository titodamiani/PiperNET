import subprocess
import argparse
from src.rnaseq_utils import get_config
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Cluster proteome using CD-HIT.")
    parser.add_argument("--config", dest='config_path', help="Path to config.yaml", default='config/rnaseq_dataprep.yaml')
    args = parser.parse_args()

    #load sample folder paths from config.yaml
    assembly_paths = get_config(args.config_path, data='assembly')

    #iterate over path dictionary
    for input in assembly_paths.values():
        
        try:
            #i/o file paths
            output = input.with_name('transcriptome_clstr.pep')

            #construct command
            min_identity = str(0.95)
            print(f'running CD-HIT on {input} with identity threshold: {min_identity}')
            command = f"cd-hit -i {input} -o {output} -c {min_identity} -n 5 -g 1 -d 0"

            #run command
            subprocess.run(command, shell=True)

        #skip to next file if input is not found
        except FileNotFoundError:
            print(f'{input} not found. Skipping to next sample.')
            continue

if __name__ == "__main__":
    main()