import argparse
from pathlib import Path
import pandas as pd
from src.lcms_utils import *



def main():
    parser = argparse.ArgumentParser(description="Merge MZmine feature table with SIRIUS predictions and GNPS results.")
    parser.add_argument("--config", dest='config_path', help="Path to config.yaml", default='config/lcms_dataprep.yaml')
    args = parser.parse_args()

    #load paths dict from config.yaml
    config_path = 'config/lcms_dataprep.yaml'
    config = get_config(config_path, data='lc-ms')


    ### Sample list ###
    metadata = pd.read_csv(config['metadata'], sep='\t')
    filenames = metadata['filename']
    sample_names = metadata['ATTRIBUTE_Replicate']


    ### MZmine ###
    ftable = pd.read_csv(config['ftable'])
    ftable.columns = [col.replace(' Peak area', '') if ' Peak area' in col else col for col in ftable.columns] #remove 'Peak area' from sample columns
    ftable = ftable[['row ID', 'row m/z', 'row retention time', 'correlation group ID', 'best ion', 'neutral M mass'] + filenames.tolist()] #keep relevant columns
    rename_map = dict(zip(filenames, sample_names)) #rename sample columns
    rename_map.update({'row ID': 'feat_ID',
                    'row m/z': 'mz',
                    'row retention time': 'RT',
                    'correlation group ID': 'corrGroup_ID',
                    'best ion':'adduct',
                    'neutral M mass': 'neutral mass'}) #rename ID columns
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
    annotations = pd.read_csv(config['annotations'], usecols=['id', 'compound_name'], na_values='')
    annotations.rename(columns={'id': 'feat_ID',
                                'compound_name': 'customDB'}, inplace=True)
    ftable = pd.merge(ftable, annotations, on='feat_ID', how='left') #create 'customDB' column
    ftable['customDB'] = ftable['customDB'].notnull() #convert to boolean


    ### SIRIUS ###
    #CSI:FingerID
    fingerid = pd.read_csv(config['fingerid'], sep='\t', usecols= ['smiles', 'name', 'CSI:FingerIDScore', 'ConfidenceScore', 'molecularFormula', 'featureId'])
    fingerid.rename(columns={'name': 'FingerID_name',
                            'ConfidenceScore':'COSMIC_score',
                            'CSI:FingerIDScore':'FingerID_score',
                            'smiles': 'FingerID_smiles',
                            'molecularFormula':'FingerID_formula',
                            'featureId':'feat_ID'}, inplace=True)

    #CANOPUS
    canopus = pd.read_csv(config['canopus'], sep='\t', usecols= ['NPC#pathway', 'NPC#pathway Probability', 'NPC#superclass', 'NPC#superclass Probability', 'NPC#class', 'NPC#class Probability', 'molecularFormula', 'featureId'])
    canopus.rename(columns={'NPC#pathway':'NPCpathway',
                            'NPC#pathway Probability':'NPCpathway_score',
                            'NPC#superclass':'NPCsuperclass',
                            'NPC#superclass Probability':'NPCsuperclass_score',
                            'NPC#class':'NPCclass',
                            'NPC#class Probability':'NPCclass_score',
                            'molecularFormula':'CANOPUSformula',
                            'featureId':'feat_ID'}, inplace=True)



    ### GNPS2 node table ###
    gnps = pd.read_csv(config['gnps2_ntable'], sep='\t')
    cols = ['cluster index', 'parent mass', 'RTMean', 'component'] + [col for col in gnps.columns if 'GNPSGROUP' in col]
    gnps = gnps[cols]
    gnps.rename(columns={'cluster index': 'feat_ID',
                    'parent mass': 'mz',
                    'RTMean': 'RT',
                    'component': 'network_ID'}, inplace=True) #rename columns

    #create 'network_size' column
    network_size = gnps.groupby('network_ID')['network_ID'].count()
    gnps = gnps.merge(network_size, how='left', left_on ='network_ID', right_index=True, suffixes = (None, '_count')).rename(columns={'network_ID_count': 'network_size'}) 


    #create 'Log2_intensity' column (using ATTRIBUTE_Replicate columns)
    rep_cols = gnps.filter(regex='ATTRIBUTE_Replicate')
    gnps['Log2-intensity'] = np.log2(rep_cols.sum(axis=1))


    #keep only 'ATTRIBUTE_Sample' columns for pie-charts mapping
    cols_to_keep = [col for col in gnps.columns if 'ATTRIBUTE' not in col or 'ATTRIBUTE_Tissue' in col]
    ntable = gnps[cols_to_keep]
    ntable.columns = ntable.columns.str.replace('ATTRIBUTE_Tissue:GNPSGROUP:', '') #shorten column names

    ### GNPS2 library matches ###
    mslib = pd.read_csv(config['gnps2_lib'], sep='\t')
    mslib = mslib[['#Scan#', 'Compound_Name', 'LibMZ', 'MassDiff', 'MQScore', 'Smiles', 'Instrument']]
    mslib.rename(columns={'#Scan#': 'feat_ID',
                    'Compound_Name': 'LibMatch_name',
                    'LibMZ': 'Lib_mz',
                    'MassDiff': 'Lib_MassDiff',
                    'MQScore': 'LibMatch_score',
                    'Instrument': 'LibMatch_instrument',
                    'Smiles': 'LibMatch_smiles'}, inplace=True) #rename columns



    ### Merge and Export ###
    #MZmine feateure table
    ftable = pd.merge(ftable, canopus, on='feat_ID', how='left')
    ftable = pd.merge(ftable, fingerid, on='feat_ID', how='left')
    ftable = pd.merge(ftable, gnps[['network_ID', 'network_size', 'feat_ID']], on='feat_ID', how='left')
    ftable.insert(6, 'network_ID', ftable.pop('network_ID'))
    ftable.insert(7, 'network_size', ftable.pop('network_size'))
    ftable.to_csv(config['ftable_merged'], index=False)

    #GNPS node table
    ntable = pd.merge(ntable, ftable[['feat_ID', 'adduct', 'neutral mass', 'customDB', 'Detected', 'corrGroup_ID']], on='feat_ID', how='left')
    ntable = pd.merge(ntable, canopus, on='feat_ID', how='left')
    ntable = pd.merge(ntable, fingerid, on='feat_ID', how='left')
    ntable = pd.merge(ntable, mslib, on='feat_ID', how='left')
    ntable['mz'] = ntable['mz'].round(4) #round m/z to 4 decimals
    ntable['RT'] = ntable['RT'].round(2) #round RT to 2 decimals
    ntable.to_csv(config['ntable_clean'], index=False)


if __name__ == "__main__":
    main()