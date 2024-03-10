import pandas as pd
import yaml
import re
import os
import urllib.request
from pathlib import Path
import subprocess
from io import StringIO
from Bio import SeqIO


# ########## Load config.yaml files ##########
# def get_config(config_path, **kwargs):
    
#     #load config.yaml as dict
#     with open(config_path, "r") as handle:
#         config = yaml.safe_load(handle)
    
#     #keep portion of config.yaml based on provided keys 
#     try:
#         for level in kwargs.values():
#             if isinstance(config[level], str):
#                 return config[level]
#             else:
#                 config = config[level]
        
#         #convert str paths to Path objects
#         config = {k: Path(v) if isinstance(v, str) else v for k, v in config.items()}
    
#         return config
#     except KeyError:
#         raise KeyError("Key(s) not found in config.yaml")

########## Load config.yaml files ##########
def get_config(config_path: str, filepaths: bool=True, **kwargs):
    """
    Load a config.yaml file as dictionary.

    Parameters:
    config_path (str): path to the YAML configuration file.
    filepaths (bool): if True, convert all string values to Path objects.
    **kwargs: keyword arguments representing keys in the YAML file.

    Returns:
    dict: config.yaml file as dictionary.

    Raises:
    KeyError: if one of the keys in kwargs is not found in the config.yaml.
    """

    #load config.yaml as dict
    with open(config_path, "r") as handle:
        config = yaml.safe_load(handle)
    
    #keep portion of config.yaml based on provided keys 
    try:
        for argument, value in kwargs.items(): #iterates over config levels provided in kwargs"
            if value not in config: #if one of the provided levels is not in the config.yaml, raise error.
                raise KeyError(f'Key {value} not found in the config.yaml') 
            if isinstance(config[value], str): #if config[value] is a string, return it and end the loop
                return config[value]
            else: #if config[value] is not a string (i.e., it's still a dict), update config to be the nested dict and go to next iteration
                config = config[value]
        
        #convert str paths to Path objects
        if filepaths:
            config = {k: Path(v) if isinstance(v, str) else v for k, v in config.items()}
    
        return config
    
    except KeyError as e:
        raise KeyError(str(e))



########## Import SonicParanoid output ##########
def import_sonicpd(file_path):
    sonicpd = pd.read_csv(file_path, sep='\t', index_col='group_id')
    sonicpd = sonicpd.applymap(lambda x: x.split(',') if isinstance(x, str) else x) #turn entries into lists (instead of strings)
    sonicpd = sonicpd.applymap(lambda x: [] if x == ['*'] else x) #replace empty entries with None instead of '*' 
    sonicpd.columns = [col.split('.')[0] for col in sonicpd.columns] #remove .pep from col names
    return sonicpd


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



########## BLASTp ##########
def run_blastp(seq_paths, blastDB_paths, min_cov = 0, min_sim = 0):
    
    '''To document'''
    
    hits = pd.DataFrame()  # empty DataFrame
    
    for query in seq_paths.values():
        for sample, blastDB in zip(blastDB_paths.keys(), blastDB_paths.values()):
            print(f'Querying {query} in {sample} transcriptome...')

            #construct command
            command = f'blastp -db {blastDB} -query {query} -evalue 1e-50 -outfmt "10 delim=, qacc sseqid qcovs score evalue pident ppos sseq"'

            #run command
            output = subprocess.run(command, shell=True, capture_output=True, text=True)

            #convert to dataframe
            if output.returncode == 0:
                output = output.stdout
                output = pd.read_csv(StringIO(output), header=None, names=['query', 'hit', 'coverage', 'score', 'evalue', 'identity', 'similarity', 'sequence'])

                #create full_id column
                output['blast_id'] = output.apply(lambda row: f"{row['hit']}_cov{row['coverage']}_sim{round(row['similarity'])}", axis=1)
                
                #print number of retrieved hits
                print(f"{output.shape[0]} hits found for {query.stem} in {sample} transcriptome.")

                #filter by coverage and similarity
                filt = (output['coverage'] > min_cov) & (output['similarity'] > min_sim)
                output = output.loc[filt]

                #append output to blast_results
                hits = pd.concat([hits, output], ignore_index=True)

            else:
                print(f"An error occurred while running the BLASTp command!")

    # #remove duplicates
    # hits['query'] = hits['query'].apply(lambda x: [x]) #convert to list
    # grouped = hits.groupby('hit').agg({'query': 'sum'})
    # hits = hits.drop(columns='query')
    # hits = pd.merge(hits, grouped, left_on='hit', right_index=True)
    # hits['query'] = hits['query'].apply(set) #convert to set

    return hits