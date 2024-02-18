import pandas as pd
import yaml
import re
import os
import urllib.request
from pathlib import Path
from Bio import SeqIO


########## Load config.yaml files ##########
def get_config(config_path, **kwargs):
    
    #load config.yaml as dict
    with open(config_path, "r") as handle:
        config = yaml.safe_load(handle)
    
    #keep portion of config.yaml based on provided keys 
    try:
        for level in kwargs.values():
            config = config[level]
        
        #convert str paths to Path objects
        config = {k: Path(v) if isinstance(v, str) else v for k, v in config.items()}

        return config
    except KeyError:
        raise KeyError("Key(s) not found in config.yaml")



########## Extract orf type from rnaSPAdes header ##########
def extract_orf_type(header):
    if 'complete' in header:
        return 'complete'
    elif '3prime_partial' in header:
        return '3prime_partial'
    elif '5prime_partial' in header:
        return '5prime_partial'
    elif 'internal' in header:
        return 'internal'
    else:
        return None
    


########## Extract orf orientation from rnaSPAdes header ##########
def extract_orientation(header):
    if '(+)' in header:
        return '(+)'
    elif '(-)' in header:
        return '(-)'
    else:
        return None



########## Extract TransDecoder score from rnaSPAdes header ##########
def extract_score(header):
    match = re.search('score=(\d+\.\d+)', header)
    if match:
        return match.group(1)
    else:
        return 0  # or any other default value you prefer



########## Shorten full rnaSPAdes ID ##########
def shorten_id(df):
    df['id'] = df['id'].str.split('_', n=6).str[6]
    return df




########## Download sequences (FASTA) in enzymeDB ##########
def download_enzymeDB(db_path, output_folder):

    enzymeDB = pd.read_csv(db_path) #import enzymeDB

    #create enzyme_seq folder if it doesn't exist
    if Path(output_folder).exists():
            print(f"'enzyme_seq' folder already exists at {output_folder}")
    else:
        output_folder.mkdir(parents=True, exist_ok=True)
        print(f"Creating 'enzyme_seq' folder at {output_folder}...")

    #download sequences (if not already present) and shorten headers
    for _, row in enzymeDB.iterrows():
        enzyme = row['Enzyme']
        file_path = os.path.join(output_folder, f'{enzyme}.fasta')

        #check if file already exists
        if os.path.exists(file_path):
            print(f"{enzyme} sequence already present at '{file_path}'")
        else:
            urllib.request.urlretrieve(url=row['URL'], filename=file_path)
            print(f"{enzyme} sequence downloaded at '{file_path}'")

        #shorten header
        for record in SeqIO.parse(file_path, 'fasta'):
            record.description = ""

            with open(file_path, 'w') as modified_fasta:
                SeqIO.write(record, modified_fasta, 'fasta')

        records = []
        with open(file_path, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                record.id = enzyme
                record.description = ''
                records.append(record)
        with open(file_path, 'w') as f:
            SeqIO.write(records, f, 'fasta')



########## Reorder tissues ##########
def reorder_tissues(original):
    tissue_order = ['leaf', 'stem', 'root', 'youfr', 'medfr', 'oldfr']
    return [tissue for tissue in tissue_order if tissue in original]