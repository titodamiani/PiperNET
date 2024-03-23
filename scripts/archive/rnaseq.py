import pandas as pd
from Bio import SeqIO
import csv
import re
import os
import urllib.request
import subprocess
from io import StringIO
import yaml



#read config.yaml
def load_config(path, level1=None, level2=None):
    with open(path, "r") as handle:
        config = yaml.safe_load(handle)
    
    if level2:
        return config[level1].get(level2, 'Level2 not found')
    else:
        return config.get(level1, 'Level1 not found1') 




########## Parse transcriptome as pandas DataFrame ##########

def parse_transcriptome(file_path):
    """Parses transXpress output files transcriptome_annotated.fasta or transcriptome.fasta (rnaSPAdes) and return pandas DataFrame:
    - Sequence ID is split into 4 columns: `protein_ID`, `Node`, `Transcript_length`, `Transcript_coverage`
    - Protein sequence is stored into the `Sequence` column
    - All the rest is store in to the `Description` column

    Parameters:
    file_path (str)

    Returns:
    transcriptome_df: pandas dataframe
   """

    #create empty dataframe
    transcriptome_df = pd.DataFrame(columns=["Long_ID", "Sequence", "Description"])
    
    #parse proteome
    with open(file_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            new_row = {"Long_ID": record.id, "Sequence": str(record.seq), "Full header": str(record.description)}
            transcriptome_df = pd.concat([transcriptome_df, pd.DataFrame([new_row])], ignore_index=True)

    #split 'ID_long' column after every '_'
    ID_splits = transcriptome_df['Long_ID'].str.split('_', expand=True) 


    #add columns to proteome_df
    transcriptome_df[['Node', 'Transcript_length', 'Transcript_coverage']]  = ID_splits[[1,3,5]]
    transcriptome_df['Gene_ID']  = ID_splits[6] + '_' + ID_splits[7]

    #reorder columns
    transcriptome_df = transcriptome_df[['Gene_ID', 'Long_ID',  'Node', 'Transcript_length', 'Transcript_coverage', 'Sequence', 'Full header']]

    return transcriptome_df





########## Parse proteome as pandas DataFrame ##########

def parse_proteome(file_path):
    """Parses transXpress output files transcriptome_annotated.pep or transcriptome.pep (rnaSPAdes) and return pandas dataframe.

    Parameters:
    file_path (str)

    Returns:
    proteome_df: pandas dataframe
   """

    #create empty dataframe
    proteome_df = pd.DataFrame(columns=["ID_long", "Sequence", "Description"])
    
    #parse proteome
    with open(file_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            new_row = {"Long_ID": record.id, "Sequence": str(record.seq), "Full header": str(record.description)}
            proteome_df = pd.concat([proteome_df, pd.DataFrame([new_row])], ignore_index=True)

    #split 'ID_long' column after every '_'
    ID_splits = proteome_df['Long_ID'].str.split('_', expand=True) 


    #add columns to proteome_df
    proteome_df[['Node', 'Transcript_length', 'Transcript_coverage']]  = ID_splits[[1,3,5]]
    proteome_df['Protein_ID']  = ID_splits[6] + '_' + ID_splits[7]

    #reorder columns
    proteome_df = proteome_df[['Protein_ID', 'Long_ID',  'Node', 'Transcript_length', 'Transcript_coverage', 'Sequence', 'Full header']]

    return proteome_df





########## Parse BLASTP results as pandas DataFrame ##########

def parse_blastp(file_path):
    """Parses transXpress output file sprotblastp_orfs.out and return a pandas Dataframe with the BLASTP hit and coresponding E-value.

    Parameters:
    file_path (str)

    Returns:
    blastp_df: pandas DataFrame
   """
    
    blastp_annotations={} #empty dict

    with open(file_path) as handle:
        csv_reader = csv.reader(handle, delimiter="\t")
        for row in csv_reader:
            if (len(row) < 13): continue
            blastp_annotations[row[0]] = [row[12], row[10]]
    

    blastp_df = pd.DataFrame.from_dict(blastp_annotations, orient='index', columns=['BLASTP hit', 'BLASTP E-value']).reset_index().rename(columns={'index': 'Long_ID'}) #convert to dataframe
    
    return blastp_df





### Parse Pfam predictions as pandas DataFrame

def parse_pfam(file_path):
    """Parses transXpress output file pfam_orfs.out and return a pandas Dataframe with the Pfam ID, domain description and coresponding E-value.

    Parameters:
    file_path (str)

    Returns:
    pfam_df: pandas DataFrame
   """
    
    pfam_annotations={} #empty dict

    with open(file_path) as handle:
        for line in handle:
            
            if (line.startswith("#")): continue
            row = re.split(" +", line, 22)
            if row[3] not in pfam_annotations:
                pfam_annotations[row[3]] = [row[1], row[22].rstrip('\n'), row[6]]
        

    pfam_df = pd.DataFrame.from_dict(pfam_annotations, orient='index', columns=['Pfam hit', 'Domain description', 'Pfam E-value']).reset_index().rename(columns={'index': 'Long_ID'}) #convert to dataframe
    
    return pfam_df





########## Merge proteome with annotations and expression table ##########

def merge_proteome_annotations(proteome, blastp, pfam, expression_table):
    """Merge preotome with, BLASTP results, Pfam predictions and expression values. All transXpress output files must be first parsed as pandas DataFrame.
        
    Parameters:
    proteome: pandas DatFrame
    blastp: pandas DatFrame
    pfam: pandas DatFrame
    expression_table: pandas DatFrame

    Returns:
    proteome: pandas DatFrame
   """
    
    #merge BLASTP results
    proteome = proteome.merge(blastp, on='Long_ID', how='left')
    
    #merge pfam results
    proteome = proteome.merge(pfam, on='Long_ID', how='left')

    #merge expression table
    proteome['temporary_id'] = proteome['Long_ID'].str.replace(r'\.[^.]+$', '', regex=True) #create temporary_id column
    proteome = proteome.merge(expression_table, left_on='temporary_id', right_index=True)
    proteome.drop(['temporary_id'], axis=1, inplace=True)    

    return proteome





########## Shorten the original rnaSPAdes headers ##########

def shorten_rnaSPAdes_headers(expression_table):
    
    """
    Shortens the original rnaSPAdes transcript IDs in the transcriptome_expression_isoform.tsv file (e.g., NODE_17296_length_3000_cov_1042.610523_g1517_i3 -> g1517_i3).

    Parameters:
    expression_table (pd.DataFrame): expression table (transcriptome_expression_isoform.tsv imported as pandas DataFrame).

    Returns:
    pd.DataFrame: expression table with shortened IDs.
    """
    
    #extract indexes and convert to pandas Series
    original_ids = pd.Series(expression_table.index)

    #shorten original IDs
    new_ids =[]

    for id in original_ids:
       
       #split every '_' and join last 2 substrings
       short_id = '_'.join(id.split('_')[-2:])

       new_ids.append(short_id)

    #update DataFrame indexes
    expression_table.index = new_ids
   
    return expression_table





########## Download FASTA files from custom database  ##########

def downloadfasta(db, filepath, foldername='sequences'):
    """
    Download FASTA files into a folder. Links are reterive from the database.

    Parameters:
    - db (pandas.DataFrame): A DataFrame containing the database with columns 'Enzyme' and 'Link'.
    - filepath (str): The path to the directory where the folder will be created.
    - foldername (str, optional): The name of the folder to store the downloaded FASTA files. 
                                  Defaults to 'sequences'.

    Returns:
    - None

    Example usage:
    downloadfasta(db, '/path/to/directory', 'fasta_files')

    """

    #create output folder
    output_dir = os.path.join(filepath, foldername)
    os.makedirs(output_dir, exist_ok=True)

    for idx, row in db.iterrows():
        filename = os.path.join(output_dir, row['Enzyme'] + '.fasta')
        urllib.request.urlretrieve(url=row['Link'], filename=filename)





########## BLASTp search ##########

def blastp(blastDB_path, query_path, min_cov = 0, min_sim = 0):
    
    '''To document'''

    #run BLASTp
    command = f'blastp -db {blastDB_path} -query {query_path} -evalue 1e-50 -outfmt "10 delim=, qacc sseqid qcovs score evalue pident ppos sseq"'
    output = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    #check if command was successful
    if output.returncode == 0:
        blast_hits = output.stdout
    else:
        print("Error:", output.stderr)
        blast_hits = None

    if blast_hits:
        blast_hits = pd.read_csv(StringIO(blast_hits), names=['query_id', 'seq_id', 'coverage', 'total score', 'e-value', 'identity', 'similarity', 'sequence'])
        
        #create full_id
        blast_hits['full_id'] = blast_hits.apply(lambda row: f"{row['seq_id']}_cov{row['coverage']}_sim{round(row['similarity'])}", axis=1)

        #reorder columns
        blast_hits = blast_hits[['query_id', 'seq_id', 'full_id', 'coverage', 'total score', 'e-value', 'identity', 'similarity', 'sequence']]

        #filter by coverage and similarity
        filt = (blast_hits['coverage'] > min_cov) & (blast_hits['similarity'] > min_sim)
        blast_hits = blast_hits.loc[filt]

    return blast_hits