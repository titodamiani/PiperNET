import pandas as pd
import re
import os
import urllib.request
from pathlib import Path
import subprocess
from io import StringIO
from Bio import SeqIO
import matplotlib.pyplot as plt



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
def run_blastp(query_paths, blastDB_paths, out_path = '', min_cov = 0, min_sim = 0, min_bitscore = 50, create_combinedID = False, log = True):
    '''
    Run BLASTp search for each query in 'query_paths' against each database in 'blastDB_paths'.
    Results are filtered based on minimum coverage, similarity, and bit score.
    Optionally creates a combined ID for each hit concatenating hit ID, coverage, and similarity.
    If 'out_path' is provided, saves results as CSV file. 


    Parameters:
    query_paths (dict): dictionary where keys are query names and values are paths to the query files.
    blastDB_paths (dict): dictionary where keys are sample names and values are paths to the BLAST databases.
    out_path (str, optional): Path where to save the hits DataFrame as a CSV file. Default is ''.
    min_cov (int, optional): Minimum query coverage to consider a hit. Default is 0.
    min_sim (int, optional): Minimum similarity to consider a hit. Default is 0.
    min_bitscore (int, optional): Minimum bit score to consider a hit. Default is 50.
    create_combinedID (bool, optional): Whether to create a combined ID for each hit. Default is False.

    Returns:
    hits (DataFrame): A DataFrame containing the hits form all BLAST searches.
    '''

    hits = pd.DataFrame()
    
    for query, query_path in query_paths.items():

            for sample, database_path in blastDB_paths.items():
                    if log:
                        print(f'Querying {query} in {sample} transcriptome...')

                    #run blastp
                    command = f'blastp -db {database_path} -query {query_path} -evalue 1e-50 -outfmt "10 delim=, qacc sseqid qcovs score evalue pident ppos sseq"'
                    out = subprocess.run(command, shell=True, capture_output=True, text=True)

                    #convert subprocess output to df
                    if out.returncode == 0:
                            out = out.stdout
                            out_df = pd.read_csv(StringIO(out), header=None, names=['query', 'hit', 'query coverage', 'bit score', 'e-value', 'identity', 'similarity', 'hit sequence'])
                            
                            #create 'sample' column
                            out_df['sample'] = sample

                            #filter by querey_coverage, similarity, bit score
                            filt = (out_df['query coverage'] > min_cov) & (out_df['similarity'] > min_sim) & (out_df['bit score'] > min_bitscore)
                            out_df = out_df.loc[filt]
                            
                            #create 'combined ID' column
                            if create_combinedID:
                                    #create combinedID column
                                    out_df['combined ID'] = out_df.apply(lambda row: f"{row['hit']}_cov{row['query coverage']}_sim{round(row['similarity'])}", axis=1)

                            #print number of retrieved hits
                            if log:
                                print(f"{out_df.shape[0]} hits for {query} in {sample} transcriptome.")

                            #append output to hits df
                            hits = pd.concat([hits, out_df], ignore_index=True)

                    else:
                            print(f"ERROR: An error occurred while running the BLASTp command!")

    #save hits to file
    if out_path:
          if os.path.exists(os.path.dirname(out_path)):
                 hits.to_csv(out_path, index=False)
          else:
            raise FileNotFoundError(f"Directory {os.path.dirname(out_path)} not found. Please provide a valid path.")

    return hits



########## Plot counts per species ##########
def plot_count_per_species(data, title='', xlabel='', ylabel='Count'):
    """
    Plots a bar chart of counts per species.

    Parameters:
    - data: pandas.Series or pandas.DataFrame with the counts data to plot.
        Expected to be a Series where index are species and values are counts.
    - title: str, plot title. Default is ''.
    - xlabel: str, label on x-axis. Default is ''.
    - ylabel: str, label on y-axis. Default is 'Count'.
    """
    plt.figure(figsize=(7,5))
    data.plot(kind='bar')
    plt.title(title, fontdict={'fontsize': 18, 'fontweight': 'bold'})
    plt.xlabel(xlabel)
    plt.ylabel(ylabel, fontdict={'fontsize': 18})
    plt.xticks(fontsize=18)
    plt.show()



########## Assign single expression value in tissues ##########
def assign_tissue_expr(proteome, tpm_cols='', tissue_order=['leaf', 'stem', 'root', 'youfr', 'medfr', 'oldfr'], metric='mean'):
    '''
    Calculate expression (max or average) of replicates per tissue.

    Parameters:
    proteome (DataFrame): input proteome dataframe. Can be proteomes_all or subsets of it. 
    tpm_cols (list, optional): list of column names containing expression data. 
                               If not provided, columns containing 'rep' in their names are used.
    tissue_order (list, optional): tissue names in the order that they should appear in the output dataframe. 
                                   Default is ['leaf', 'stem', 'root', 'youfr', 'medfr', 'oldfr'].
    metric (str): metric used for calculating expression value across tissues. 
                  Either 'max' or 'mean'. Default is 'mean'.

    Returns:
    DataFrame: df with same index as input data and expression (max or average) for each tissue, 
               with tissues ordered according to tissue_order.
    
    Raises:
    ValueError: If an invalid metric is passed. Only 'max' and 'mean' are valid.
    '''
    
    if not tpm_cols:
        tpm_cols=[col for col in proteome.columns if 'rep' in col]

    expr_data = proteome[tpm_cols] #keep columns with expression data
    tissue_names = expr_data.columns.str.split('_').str[0] #extract tissues from column names

    if metric == 'max':
        expr_data_result = expr_data.groupby(by=tissue_names, axis=1).max().round(1) #group-by tissue, calculate max

    elif metric == 'mean':
        expr_data_result = expr_data.groupby(by=tissue_names, axis=1).mean().round(1) #group-by tissue, calculate mean

    else:
        raise ValueError("Invalid metric. Expected 'max' or 'mean'.")

    expr_data_result = expr_data_result.reindex(columns=tissue_order) #reorder columns
    return expr_data_result



########## Assign single expression value in orthogroups table ##########
def assign_orth_expr(orthologs, proteomes, tpm_cols='', metric='max'):
    '''
    Assign expression values to each entry in an orthogroups table. Expression values are retrieved from the proteomes_all df. Either maximum or average expression of orthologs across tissus can be used ('metric').
    
    Parameters:
    orthologs (list): list of orthologs (i.e., each entry in orthogroups table).
    proteomes (pd.DataFrame): merged proteomes dataframe.
    tpm_cols (str): list of columns containing expression data. Default is empty. If empty, all columns containing 'rep' are used.
    metric (str): metric used for assigning a single expression value across tissues. 
                  Either 'max' or 'mean'. Default is 'max'.
    
    Returns:
    expr (float): single expression value.
    
    Raises:
    ValueError: If an invalid metric is passed. Only 'max' and 'mean' are valid.
    '''

    if not tpm_cols:
        tpm_cols=[col for col in proteomes.columns if 'rep' in col]

    if metric == 'max':
        expr = proteomes.loc[:, tpm_cols].loc[orthologs].max().max()

    elif metric == 'mean':
        expr = proteomes.loc[:, tpm_cols].loc[orthologs].mean().mean()
        expr = round(expr, 1)

    else:
        raise ValueError("Invalid metric. Expected 'max' or 'mean'.")
    return expr



########## Assign single length value in orthogroups table ##########
def assign_len_orth(orthologs, proteomes, length_col='length', metric='max'):

    '''
    Assign lenght to each entry in orthogroups table. Sequence lenghts are retrieved from the proteomes_all df. Either maximum or average length of orthologs can be used ('metric').
    
    Parameters:
    orthologs (list): list of orthologs (i.e., each entry in orthogroups table).
    proteomes (pd.DataFrame): merged proteomes dataframe.
    length_col (str): name of column containing length data. Default is 'length'.
    metric (str): metric used for assigning single legth value. 
                  Either 'max' or 'mean'. Default is 'max'.
    
    Returns:
    length (int): single length value.
    
    Raises:
    ValueError: If an invalid metric is passed. Only 'max' and 'mean' are valid.
    '''

    if metric == 'max':
        length = proteomes[length_col].loc[orthologs].max()
    elif metric == 'mean':
        length = proteomes[length_col].loc[orthologs].mean()
    else:
        raise ValueError("Invalid metric. Expected 'max' or 'mean'.")
    return length



########## Plot ORF type and sequence lenght distributions ##########
def plot_orftype_and_lenght_distribution(data, colormap=None, bin_size=25):
    """
    Generate two plots: 1) pie chart of TransDecoder ORF type count; 2) Sequence lenght distribution as stacked bar plot

    Parameters:
    data (DataFrame): proteomes_all dataframe, or subsets of it.
    colormap (dict, optional): A dictionary mapping ORF types to colors. If None, a default colormap is used.
    bin_size (int, optional): The size of the bins for the sequence length distribution. Default is 25.

    Returns:
    None. The function directly plots the distributions using matplotlib.
    """
    
    #count orf_type
    orf_type_count = data['orf_type'].value_counts()

    #sequence length distribution
    bins = pd.cut(data['length'], bins=range(0, data['length'].max(), bin_size)) #bin data['length']
    count = data.groupby([bins, 'orf_type']).size().unstack() #count occurrence of orf_types in each bin
    count = count[orf_type_count.index.tolist()] #reorder according to orf_type_count
    count = count[count.sum(axis=1) > 0] #drop all-zeros rows

    #color mapping
    if colormap is None:
        colormap = {'complete': 'green',
                    '5prime_partial': 'orange',
                    '3prime_partial': 'lightcoral',
                    'internal': 'red'}
    colors=orf_type_count.index.map(lambda x: colormap[x])

    #create figures
    fig, axs = plt.subplots(1, 2, figsize=(15, 7))

    #pie chart
    orf_type_count.plot(kind='pie', autopct='%1.1f%%', colors=colors, legend=False, labels=None, textprops={'fontsize': 18}, ax=axs[0])
    axs[0].set_title('ORF type', fontdict={'fontsize': 18, 'fontweight': 'bold'})
    axs[0].set_ylabel('')  # remove y-label
    axs[0].legend(orf_type_count.index, bbox_to_anchor=(1, 1), loc='upper left') #legend
    
    #stacked bar plot
    count.plot(kind='bar', stacked=True, color=colors, ax=axs[1]) # color=colors,
    axs[1].set_title('Sequence lenght distribution', fontdict={'fontsize': 18, 'fontweight': 'bold'})
    axs[1].tick_params(axis='x', labelsize=14)  # increase x-tick labels font size
    axs[1].legend(title=None, loc='upper left')  # remove legend title

    plt.tight_layout()
    plt.show()