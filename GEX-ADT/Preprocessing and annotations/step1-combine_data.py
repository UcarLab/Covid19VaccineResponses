import muon as mu
import scanpy as sc
import seaborn as sns
from matplotlib import pyplot as plt
from venn import venn

import pandas as pd
import numpy as np
import os
from os import path
from glob import glob

import anndata


set1_pool1 = {
    'hash11' : '1002-V1',
    'hash12' : '1002-V2',
    'hash14' : '1002-V3',
    'hash15' : '1002-V4',
    'hash16' : '1002-V5',
    'hash18' : '1002-V6',
    'hash19' : '1002-V7',
    'hash21' : '1002-V8',
}


set1_pool2 = {
    'hash11' : '1022-V1',
    'hash12' : '1022-V2',
    'hash14' : '1022-V3',
    'hash15' : '1022-V4',
    'hash16' : '1022-V5',
    'hash18' : '1022-V6',
    'hash19' : '1022-V7',
    'hash21' : '1022-V8',
}


set1_pool3 = {
    'hash11' : '2005-V1',
    'hash12' : '2005-V2',
    'hash14' : '2005-V3',
    'hash15' : '2005-V6',
    'hash16' : '2005-V7',
    'hash18' : '2005-V8',
}



set2_pool1 = {
    'hash11' : '1005-V1',
    'hash12' : '1005-V2',
    'hash13' : '1005-V4',
    'hash14' : '1005-V5',
    'hash15' : '1005-V6',
    
    'hash16' : '1010-V1',
    'hash17' : '1010-V2',
    'hash18' : '1010-V4',
    'hash19' : '1010-V5',
    'hash20' : '1010-V6',
}


set2_pool2 = {
    'hash11' : '1024-V1',
    'hash12' : '1024-V2',
    'hash13' : '1024-V4',
    'hash14' : '1024-V5',
    'hash15' : '1024-V6',
    
    'hash16' : '1015-V1',
    'hash17' : '1015-V2',
    'hash18' : '1015-V4',
    'hash19' : '1015-V5',
    'hash20' : '1015-V6',
}


set2_pool3 = {
    'hash11' : '1028-V1',
    'hash12' : '1028-V2',
    'hash13' : '1028-V4',
    'hash14' : '1028-V5',
    'hash15' : '1028-V6',
    
    'hash16' : '1021-V1',
    'hash17' : '1021-V2',
    'hash18' : '1021-V4',
    'hash19' : '1021-V5',
    'hash20' : '1021-V6',
}


set2_pool4 = {
    'hash11' : '1031-V1',
    'hash12' : '1031-V2',
    'hash13' : '1031-V4',
    'hash14' : '1031-V5',
    'hash15' : '1031-V6',
    
    'hash16' : '2004-V1',
    'hash17' : '2004-V2',
    'hash18' : '2004-V4',
    'hash19' : '2004-V5',
    'hash20' : '2004-V6',
}


set2_pool5 = {
    'hash11' : '1011-V1',
    'hash12' : '1011-V2',
    'hash13' : '1011-V3',
    
    'hash14' : '2007-V1',
    'hash15' : '2007-V2',
    'hash16' : '2007-V3',
}


set2_pool6 = {
    'hash11' : '2001-V1',
    'hash12' : '2001-V2',
    'hash13' : '2001-V4',
    'hash14' : '2001-V5',
    'hash15' : '2001-V6',
    
    'hash16' : '2008-V1',
    'hash17' : '2008-V2',
    'hash21' : '2008-V4',
    'hash19' : '2008-V5',
    'hash20' : '2008-V6',
}


run_hashing_schema = {
    '2022-03-16-1-1' : set1_pool1,
    '2022-03-16-1-2' : set1_pool1,
    '2022-03-16-1-3' : set1_pool1,
    '2022-03-16-1-4' : set1_pool1,
    
    
    '2022-03-16-1-5' : set1_pool2,
    '2022-03-16-1-6' : set1_pool2,
    '2022-03-16-1-7' : set1_pool2,
    '2022-03-16-1-8' : set1_pool2,
    
    
    '2022-03-16-2-1' : set1_pool3,
    '2022-03-16-2-2' : set1_pool3,
    '2022-03-16-2-3' : set1_pool3,


    '2023-02-20-1-1' : set2_pool1,
    '2023-02-20-1-2' : set2_pool1,
    '2023-02-20-1-3' : set2_pool1,

    '2023-02-20-1-5' : set2_pool2,
    '2023-02-20-1-6' : set2_pool2,
    '2023-02-20-1-7' : set2_pool2,
    
    '2023-02-21-1-2' : set2_pool3,
    '2023-02-21-1-3' : set2_pool3,
    '2023-02-21-1-4' : set2_pool3,
    
    
    '2023-02-21-1-5' : set2_pool4,
    '2023-02-21-1-6' : set2_pool4,
    '2023-02-21-1-7' : set2_pool4,
    
    
    '2023-03-01-1-3' : set2_pool5,
    '2023-03-01-1-4' : set2_pool5,
    '2023-03-01-2-2' : set2_pool5,
    
    
    '2023-03-02-1-1' : set2_pool6,
    '2023-03-02-1-2' : set2_pool6,
    '2023-03-02-1-3' : set2_pool6
}




run_to_pool = {
    '2022-03-16-1-1' : 'set1_pool1',
    '2022-03-16-1-2' : 'set1_pool1',
    '2022-03-16-1-3' : 'set1_pool1',
    '2022-03-16-1-4' : 'set1_pool1',
    
    
    '2022-03-16-1-5' : 'set1_pool2',
    '2022-03-16-1-6' : 'set1_pool2',
    '2022-03-16-1-7' : 'set1_pool2',
    '2022-03-16-1-8' : 'set1_pool2',
    
    
    '2022-03-16-2-1' : 'set1_pool3',
    '2022-03-16-2-2' : 'set1_pool3',
    '2022-03-16-2-3' : 'set1_pool3',


    '2023-02-20-1-1' : 'set2_pool1',
    '2023-02-20-1-2' : 'set2_pool1',
    '2023-02-20-1-3' : 'set2_pool1',

    '2023-02-20-1-5' : 'set2_pool2',
    '2023-02-20-1-6' : 'set2_pool2',
    '2023-02-20-1-7' : 'set2_pool2',
    
    '2023-02-21-1-2' : 'set2_pool3',
    '2023-02-21-1-3' : 'set2_pool3',
    '2023-02-21-1-4' : 'set2_pool3',
    
    
    '2023-02-21-1-5' : 'set2_pool4',
    '2023-02-21-1-6' : 'set2_pool4',
    '2023-02-21-1-7' : 'set2_pool4',
    
    
    '2023-03-01-1-3' : 'set2_pool5',
    '2023-03-01-1-4' : 'set2_pool5',
    '2023-03-01-2-2' : 'set2_pool5',
    
    
    '2023-03-02-1-1' : 'set2_pool6',
    '2023-03-02-1-2' : 'set2_pool6',
    '2023-03-02-1-3' : 'set2_pool6',
}

multiome = mu.read_10x_h5('aggregated_samples/outs/filtered_feature_bc_matrix.h5')

htos = glob('KITE/HTO/Demuxed Results/*-demuxed_hashes.csv')
htos = {path.basename(p)[:len('xxxx-xx-xx-x-x')]: pd.read_csv(p, index_col=0, header = 0, names = ['hash']) for p in htos}


sample_number_to_id = pd.read_csv('aggregation.csv', usecols=['library_id']).iloc[:, 0].to_dict()
sample_number_to_id = {key+1:value for key, value in sample_number_to_id.items()}
id_to_sample_number = {value:key for key, value in sample_number_to_id.items()}

for run, df in htos.items():
    mapping = run_hashing_schema[run]
    df['person_visit'] = df['hash'].map(lambda hash: mapping.get(hash, hash))

multiome.obs['run_number'] = multiome.obs_names.str.slice(start=17).astype(int)

multiome.obs['library'] = multiome.obs['run_number'].apply(lambda x: sample_number_to_id[x]).values


for run, htos_df in htos.items():
    sample_number = id_to_sample_number[run]
    htos_df['aggregated_index'] = htos_df.index.str.replace('1', str(sample_number))


multiome.obs['person_visit'] = 'Unknown'
for run, htos_df in htos.items():
    multiome.obs.loc[htos_df.aggregated_index.values, 'person_visit'] = htos_df.person_visit.values

def extract_person(person_visit):
    if '-' in person_visit:
        return person_visit[:4]
    else:
        return person_visit

multiome.obs['person'] = multiome.obs.person_visit.apply(extract_person)

def extract_visit(person_visit):
    if '-' in person_visit:
        return person_visit[-2:]
    else:
        return person_visit

multiome.obs['visit'] = multiome.obs.person_visit.apply(extract_visit)


print(multiome.obs.head())


#ADT

def kite_output_as_adata(folder):
    adt = sc.read_mtx(f'{folder}/matrix.mtx.gz')
    adt = adt.transpose()
    features = pd.read_csv(f'{folder}/features.tsv.gz', header=None)[0]
    barcodes = pd.read_csv(f'{folder}/barcodes.tsv.gz', header=None)[0]
    adt.var_names = features
    adt.obs_names = barcodes
    return adt


adt_folders = glob('KITE/ADT/2023*') + glob('KITE/ADT/2022*')

adts = {}
for adt_folder in adt_folders:
    run = path.basename(adt_folder)
    hto_data = htos[run]
    hto_valid_i = hto_data[hto_data.hash.str.startswith('hash')].index
    adt_data = kite_output_as_adata(adt_folder)
    adt_data.obs_names = adt_data.obs_names.to_series() + '-1'
    adt = adt_data[hto_valid_i].copy()
    adts[run] = adt
    
total = 0
for run, adt in adts.items():
    sample_number = id_to_sample_number[run]
    print(run, adt.shape)
    total+=adt.shape[0]
print(total)


# adjusting barcode suffixes of adts
for run, adt in adts.items():
    sample_number = id_to_sample_number[run]
    adt.obs_names = adt.obs_names.str.replace('1', str(sample_number))


adt_merged = anndata.concat(adts, join='outer', fill_value=0)

# Note!:
# 'AC-CD158a-A0420' and 'AC-CD158-A0420' are two different names for THE SAME THING! Unfortunately different names used in the two sets (the initial and the second)
adt_merged.X[:, adt_merged.var_names=='AC_CD158_A0420'] += adt_merged.X[:, adt_merged.var_names=='AC_CD158a_A0420']
adt_merged = adt_merged[:, ~(adt_merged.var_names=='AC_CD158a_A0420')].copy()

# In the pilot study they used the antibody AC_CD14_A0051 for CD14 and it actually perfromed so poorly.
# They replaced it with AC_CD141_A0163 which worked much better. We add the counts of the first to the latter and remove the first
adt_merged.X[:, adt_merged.var_names=='AC_CD141_A0163'] += adt_merged.X[:, adt_merged.var_names=='AC_CD14_A0051']
adt_merged = adt_merged[:, ~(adt_merged.var_names=='AC_CD14_A0051')].copy()


multiome_adt_common_index = multiome.obs_names.intersection(adt_merged.obs_names)

multiome = multiome[multiome_adt_common_index].copy()

mudata = multiome.copy()
mudata.mod['adt'] = adt_merged.copy()
mudata = mu.MuData(mudata.mod)

mudata.obs[multiome.obs.columns] = multiome.obs.copy()


controls_dict = {
    'AC_IgG1_A0090' : 'Cntrl_IgG1_1',
    'AC_IgG1_A0236' : 'Cntrl_IgG1_2',
    'AC_IgG2a_A0091' : 'Cntrl_IgG2a_1', 
    'AC_IgG2a_A0238' : 'Cntrl_IgG2a_2', 
    'AC_IgG2b_A0092' : 'Cntrl_IgG2b_1', 
    'AC_IgG2b_A0095' : 'Cntrl_IgG2b_2'
}
new_var_names = [f'{name[3:-6]}_prot' if "IgG" not in name else controls_dict[name] for name in mudata['adt'].var_names]
mudata['adt'].var_names = new_var_names
mudata.update()


not_incommon = {'CD14_prot', 'CD112_prot', 'CD155_prot', 'CD244_prot', 'CD7_prot', 'IgM_prot'}
controls = {'Cntrl_IgG1_1', 'Cntrl_IgG1_2', 'Cntrl_IgG2a_1', 'Cntrl_IgG2a_2', 'Cntrl_IgG2b_1', 'Cntrl_IgG2b_2'}
mudata['adt'].var['is_in_common'] = ~mudata['adt'].var_names.isin(not_incommon)
mudata['adt'].var['is_control'] = mudata['adt'].var_names.isin(controls)

mudata.update()

mudata.var_names_make_unique()

#Metadata

metadata = pd.read_csv('Antibody titers and metadata.csv')
metadata.head()
metadata = metadata.fillna(0)
metadata = metadata.rename(columns={'VACCINE':'Vaccine', 'AGE (yrs)' : 'Age', })
metadata['baseline_titer']=  metadata['V1 (Prevaccination)']
metadata['final_titer'] = metadata['V7 (2nd vaccine Day 70)'] + metadata['V6 (1st vaccine Day 35)']
metadata['fold_change'] = metadata['final_titer']/metadata['baseline_titer']
metadata['SUBJECT ID'] = metadata['SUBJECT ID'].astype(str) 
metadata.drop(columns=['V1 (Prevaccination)', 'V4 (1st vaccine Day 25)', 
                       'V6 (1st vaccine Day 35)', 'V7 (2nd vaccine Day 70)'], inplace=True)

mudata.obs = mudata.obs.reset_index().merge(metadata, left_on = 'person', right_on = 'SUBJECT ID', how='left').set_index('index')

mudata.obs['Visit'] = mudata.obs['visit'].map({
    'V1': 'Baseline 1',
    'V2': 'Vac1 D1',
    'V3': 'Vac1 D7',
    'V4': 'Baseline 2',
    'V5': 'Vac2 D1',
    'V6': 'Vac2 D7',
    'V7': 'Day70', # From the start
    'V8': 'Day180', # From the start
})

mudata.obs['pool'] = mudata.obs['library'].map(run_to_pool)

mudata.obs['set'] = mudata.obs.library.str.startswith('2022').map({True:'Initial Set', False:'Second Set'})

mudata.obs['barcode'] = mudata.obs_names.copy()

print(mudata.obs.head())


mudata.write_h5mu('combined_raw_mudata.h5mu')



