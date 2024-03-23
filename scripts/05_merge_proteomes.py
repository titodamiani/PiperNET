import argparse
import pandas as pd
from src.rnaseq_utils import get_config

def main():
    parser = argparse.ArgumentParser(description="Merge individual proteomes_csv into one and store as CSV file.")
    parser.add_argument("--config", dest='config_path', help="Path to config.yaml", default='config/rnaseq_dataprep.yaml')
    parser.add_argument("--output", dest='out_path', help="Outpath for merged proteomes_all.csv", default='data/rna-seq/proteomes_all.csv')
    args = parser.parse_args()


    #load proteomes_csv paths from config.yaml
    proteome_list = get_config(args.config_path, data='proteome_csv')
    prefixes = get_config(args.config_path, data='prefix', filepaths=False)

    # Create a new dictionary to hold the dataframes
    proteome_dfs = {}

    for sample, path in proteome_list.items(): 
        proteome_dfs[sample] = pd.read_csv(path, index_col='id')

    # Create a unique df with all proteomes
    proteomes_all = pd.DataFrame()

    for sample, proteome in proteome_dfs.items():
        proteome['sample'] = sample
        proteome.columns = [col.split('_', 1)[1] if 'rep' in col else col for col in proteome.columns]
        proteomes_all = pd.concat([proteomes_all, proteome])

    # Replace NaN with 0
    tpm_cols = [col for col in proteomes_all.columns if 'rep' in col]
    proteomes_all[tpm_cols] = proteomes_all[tpm_cols].fillna(0)

    #save to csv
    proteomes_all.to_csv(args.out_path, index=True, index_label='id')

if __name__ == "__main__":
    main()