import argparse
import yaml
from pathlib import Path
import pandas as pd
import numpy as np
import networkx as nx

config_path = Path('../config.yaml')

#read config.yaml
def load_config(section='lc-ms'):
    with open(config_path, "r") as handle:
        config = yaml.safe_load(handle)
    return config.get(section, {})  

config = load_config()

##### GNPS node table #####

#import graphml
graphml_path = Path(config.get('gnps_graphml', ''))
graph = nx.read_graphml(graphml_path)
node_table = pd.DataFrame.from_dict({node_id: attributes for node_id, attributes in graph.nodes(data=True)}, orient='index')



# #keep relevant columns
# cols = ['GNPSGROUP:' + sample for sample in sample_names] + ['cluster index', 'componentindex', 'precursor mass', 'parent mass', 'sum(precursor intensity)', 'Compound_Name', 'SharedPeaks', 'Adduct', 'neutral M mass', 'Smiles', 'Best Ion', 'IIN Best Ion=Library Adduct', 'MZErrorPPM', 'MQScore', 'MassDiff', 'RTMean'] + ['Analog:Compound_Name', 'Analog:Smiles', 'Analog:SharedPeaks', 'Analog:MassDiff', 'Analog:Adduct']
# node_table = node_table[cols]

# #rename columns
# node_table.columns = [col.replace('GNPSGROUP:', '') for col in node_table.columns]
# rename_map= {'cluster index': 'feat_ID',
#     'componentindex': 'network_ID',
#     'RTMean': 'RT',
#     'sum(precursor intensity)': 'intensity',
#     'Compound_Name': 'name_gnps',
#     'Adduct': 'adduct_gnps',
#     'Smiles': 'smiles_gnps',
#     'neutral M mass': 'neutral mass',
#     'Best Ion': 'adduct_iin',
#     'IIN Best Ion=Library Adduct': 'addcut_iin=adduct_gnps',
#     'MZErrorPPM': 'errorPPM_gnps',
#     'MassDiff': 'errorMZ_gnps'}

# node_table.rename(columns=rename_map, inplace=True)


# #create log2- and log10-intensity columns
# node_table['Log2_intensity'] = np.log2(node_table['intensity'])
# node_table['Log10_intensity'] = np.log10(node_table['intensity'])
# node_table

# #create network_zie column
# network_size = node_table.groupby('network_ID')['network_ID'].count()
# node_table = node_table.merge(network_size, how='left', left_on ='network_ID', right_index=True, suffixes = (None, '_count'))
# node_table.rename(columns={'network_ID_count': 'network_size'}, inplace=True)




# TODO: Fix graphml update (Roman)

# for node, data in graph.nodes(data=True):

#     new_data = node_table.loc[node_table['feat_ID'] == int(node)].iloc[0].to_dict()
#     nx.set_node_attributes(graph, {node: new_data})


# print(graph.nodes(data=True))

# nx.write_graphml(graph, '../20230929_piper_screen_final/data/test.graphml')



# TODO: Create `Unique network` column (Roman) and merge to MZmine feature table and GNPS node table

# def agg_samples(samples_in_network, uniqueness=0.8):
    
#     counts = Counter([len(list(s)) for s in list(samples_in_network)])
    
#     return counts[1] / sum(counts.values()) > uniqueness

# unique_networks = ftable.groupby(by='network_ID').agg({'Detected': agg_samples}).rename(columns={'Detected': 'Unique network'})


# #merge to MZmine feature_table and change column order
# ftable = ftable.merge(unique_networks, on='network_ID')
# ftable.insert(8, 'Unique network', ftable.pop('Unique network'))


# #merge to GNPS node_table
# fbmn = fbmn.merge(unique_networks, how='left', left_on ='componentindex', right_index=True) #create network_size column in the node_table