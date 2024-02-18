import argparse
from pathlib import Path
import pandas as pd
from src.lcms import *



def main():
    parser = argparse.ArgumentParser(description="Merge MZmine feature table with SIRIUS predictions and GNPS results.")
    parser.add_argument("--config", dest='config_path', help="Path to config.yaml", default='config.yaml')
    args = parser.parse_args()

    #load paths dict from config.yaml
    config = load_config(args.config_path, level1='lc-ms')

    ### Sample list ###
    metadata_path = Path(config.get('metadata', ''))
    metadata = pd.read_csv(metadata_path, sep='\t')
    filenames = metadata['filename']
    sample_names = metadata['ATTRIBUTE_Replicate']


    ### MZmine ###
    ftable_path = Path(config.get('ftable', ''))
    ftable = pd.read_csv(ftable_path)
    ftable.columns = [col.replace(' Peak area', '') if ' Peak area' in col else col for col in ftable.columns] #remove 'Peak area' from sample columns
    ftable = ftable[['row ID', 'row m/z', 'row retention time', 'correlation group ID', 'best ion'] + filenames.tolist()] #keep relevant columns
    rename_map = dict(zip(filenames, sample_names)) #rename sample columns
    rename_map.update({'row ID': 'feat_ID',
                    'row m/z': 'mz',
                    'row retention time': 'RT',
                    'correlation group ID': 'corrGroup_ID',
                    'best ion':'adduct'}) #rename ID columns
    ftable.rename(columns=rename_map, inplace=True)

    #create `corrGroup_size` column
    corrGroup_size = ftable.groupby('corrGroup_ID')['corrGroup_ID'].count() #count features in corrGroups
    ftable = ftable.merge(corrGroup_size, how='left', left_on ='corrGroup_ID', right_index=True, suffixes = (None, '_size'))
    ftable.insert(4, 'corrGroup_size', ftable.pop('corrGroup_ID_size')) #change column order

    #create 'Detected' column
    def where_detected(feature_table, sample_names):
        detected = set()
        for feat_id, intensity in enumerate(feature_table[sample_names]):
            if intensity > 0:
                detected.add(sample_names[feat_id])
        return detected

    ftable['Detected'] = ftable.apply(lambda row: where_detected(row, sample_names), axis=1)


    #import localDB annotations
    annotations_path = Path(config.get('annotations', ''))
    annotations = pd.read_csv(annotations_path, usecols=['id', 'compound_name'], na_values='')
    annotations.rename(columns={'id': 'feat_ID',
                                'compound_name': 'customDB'}, inplace=True)
    ftable = pd.merge(ftable, annotations, on='feat_ID', how='left') #create 'customDB' column
    ftable['customDB'] = ftable['customDB'].notnull() #convert to boolean


    ### SIRIUS ###
    #CSI:FingerID
    fingerid_path = Path(config.get('fingerid', ''))
    fingerid = pd.read_csv(fingerid_path, sep='\t', usecols= ['smiles', 'name', 'CSI:FingerIDScore', 'ConfidenceScore', 'molecularFormula', 'featureId'])
    fingerid.rename(columns={'name': 'FingerID_name',
                            'ConfidenceScore':'COSMIC_score',
                            'CSI:FingerIDScore':'FingerID_score',
                            'smiles': 'FingerID_smiles',
                            'molecularFormula':'FingerID_formula',
                            'featureId':'feat_ID'}, inplace=True)

    #CANOPUS
    canopus_path = Path(config.get('canopus', ''))
    canopus = pd.read_csv(canopus_path, sep='\t', usecols= ['NPC#pathway', 'NPC#pathway Probability', 'NPC#superclass', 'NPC#superclass Probability', 'NPC#class', 'NPC#class Probability', 'molecularFormula', 'featureId'])
    canopus.rename(columns={'NPC#pathway':'NPCpathway',
                            'NPC#pathway Probability':'NPCpathway_score',
                            'NPC#superclass':'NPCsuperclass',
                            'NPC#superclass Probability':'NPCsuperclass_score',
                            'NPC#class':'NPCclass',
                            'NPC#class Probability':'NPCclass_score',
                            'molecularFormula':'CANOPUSformula',
                            'featureId':'feat_ID'}, inplace=True)


    ### GNPS node table ###
    gnps_path = Path(config.get('gnps_ntable', ''))
    gnps = pd.read_csv(gnps_path, sep='\t')
    cols = ['GNPSGROUP:' + sample for sample in sample_names] + ['cluster index', 'componentindex','sum(precursor intensity)']
    gnps = gnps[cols]
    gnps.rename(columns={'cluster index': 'feat_ID',
                                'componentindex': 'network_ID',
                                'sum(precursor intensity)': 'sum_intensity'}, inplace=True)


    #create 'network_size' column
    network_size = gnps.groupby('network_ID')['network_ID'].count()
    gnps = gnps.merge(network_size, how='left', left_on ='network_ID', right_index=True, suffixes = (None, '_count'))
    gnps.rename(columns={'network_ID_count': 'network_size'}, inplace=True)


    ### Merge and Export ###
    merged_ftable_path = Path(config.get('merged_ftable', ''))
    ftable = pd.merge(ftable, canopus, on='feat_ID', how='left')
    ftable = pd.merge(ftable, fingerid, on='feat_ID', how='left')
    ftable = pd.merge(ftable, gnps[['network_ID', 'network_size', 'feat_ID']], on='feat_ID', how='left')
    ftable.insert(6, 'network_ID', ftable.pop('network_ID'))
    ftable.insert(7, 'network_size', ftable.pop('network_size'))
    ftable.to_csv(merged_ftable_path, index=False)


if __name__ == "__main__":
    main()