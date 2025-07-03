# Loading libraries

import pandas as pd
import numpy as np

# Setting file paths 

operons_file = "/space24/PGC/emhernan/3_TRN/operon_mapper/output/list_of_operons_546291"
trn_first_gene_file = "/space24/PGC/emhernan/3_TRN/output/oTF-TGfirstGene-TRN-bgM1.tsv"
orthologous_gene_file = "/space24/PGC/emhernan/3_TRN/OrthologousTFInfo/Orthologous_table.tsv"
network_info_file = "/space24/PGC/emhernan/3_TRN/networkInfo/network_tf_gene.txt"

# Reading files

# 1)  Reading and formating operon file

operons = pd.read_csv(operons_file, sep="\t")
operons["Operon"] = operons["Operon"].ffill().astype(int)
operons = operons.loc[~operons['IdGene'].isna()]


# 2) Reading and formating trn_firstGene file

trn_first_gene = pd.read_csv(trn_first_gene_file, sep="\t")
trn_first_gene = trn_first_gene[['TF_locusTag','TFname', 'TF_oldRZlocusTag', 'TG_oldRZlocusTag']]
trn_first_gene = trn_first_gene.rename(columns={'TF_locusTag': 'TFlocusTag', 'TF_oldRZlocusTag': 'TFoldRZlocusTag','TG_oldRZlocusTag': 'firstGene'})

# 3) Reading and formating orthologous genes file

orthologous_gene = pd.read_csv(orthologous_gene_file, sep="\t")
orthologous_gene = orthologous_gene[['EC_locusTag', 'NCBI_name', 'Regulondb_name', 'Abasy_name', 'Ecocyc_name', 'Synonyms', 'RegulondbID', 'RZ_locusTag_old', 'RZ_locusTag']]
orthologous_gene = orthologous_gene.drop_duplicates(subset=None, keep='first', inplace=False)

# Reading network info

network_info = pd.read_csv(network_info_file, sep="\t", comment = "#", header = None)
network_info = network_info[[0,1,2,3]]

## Setting all my functions

def get_first_gene(group):
    strand = group['Strand'].iloc[0]
    if strand == '+':
        idx = group['PosLeft'].idxmin()
    else:
        idx = group['PosLeft'].idxmax()
    return group.loc[idx, 'IdGene']


def sort_operon_genes(group):
    strand = group['Strand'].iloc[0]
    return group.sort_values(by='PosLeft', ascending=(strand == '+'))

def summarize(group):
    sorted_group = sort_operon_genes(group)
    return pd.Series({
        'TUid': group.name,
        'firstGene': get_first_gene(sorted_group),
        'TGs': ','.join(sorted_group['IdGene']),
        'No.TGs': len(sorted_group),
        'Strand': sorted_group['Strand'].unique()[0]
    })

def build_operon_summary(df):
    summary = df.groupby('Operon').apply(summarize).reset_index(drop=True)
    return summary


def summarize_tf (group):
    orto = ",".join(group['isTgOrthologous']).count('T')
    return pd.Series({
    'TFname' : group.name,
    'number_of_tgs': sum(group['No.TGs']),
    'number_orthologous_tgs': orto
    })


def setting_orthologous_tgs(data, orto_gene):
    orto_old_set = set(orto_gene['RZ_locusTag_old'])
    orto_new_set = set(orto_gene['RZ_locusTag'])
    is_orto = []
    
    for locus in map(str.strip, str(data).split(',')):
        if locus in orto_old_set or locus in orto_new_set:
            is_orto.append('T')
        else:
            is_orto.append('F')
    return pd.Series([','.join(is_orto)])

def setting_orthologous_unqiue_tgs(tgs, orto_gene):
    orto_old_set = set(orto_gene['RZ_locusTag_old'])
    orto_new_set = set(orto_gene['RZ_locusTag'])
    is_orto = []
    
    for locus in map(str.strip, str(tgs).split(',')):
        if locus in orto_old_set or locus in orto_new_set:
            is_orto.append(True)
    return len(is_orto)

def get_orthologous_interactions_and_log(row, orthologous_gene, network_info, tf_column='TFregulonDBid'):
    tf_regulonId = row[tf_column]
    tf_name = row.get('TFname', tf_regulonId)  # Use TFname if exists, and if not the ID in RegulonDB
    isOrthologousInteraction = []
    orthologous_hits = []  # Here we save the valid orthologous interactions

    for tg in row['TGs'].split(','):
        found_match = False

        for col in ['RZ_locusTag_old', 'RZ_locusTag']:
            matches = orthologous_gene[orthologous_gene[col] == tg]['RegulondbID'].values
            if len(matches):
                target_id = matches[0]
                interaction = network_info[(network_info[0] == tf_regulonId) & (network_info[2] == target_id)]
                if not interaction.empty:
                    isOrthologousInteraction.append('T')
                    orthologous_hits.append({'TFname': tf_name, 'TG_locusTag': tg})
                    found_match = True
                    break

        if not found_match:
            isOrthologousInteraction.append('F')

    return pd.Series([','.join(isOrthologousInteraction)]), orthologous_hits

def process_and_log(row):
    interactions, hits = get_orthologous_interactions_and_log(row, orthologous_gene, network_info)
    all_orthologous_hits.extend(hits)
    return interactions

########################################################################
############## Merging trn_first_gene, tu_firstGene_tgs ################
########################################################################

tu_firstGene_tgs = build_operon_summary(operons)
full_annotation = pd.merge(trn_first_gene, tu_firstGene_tgs, how = 'left', on = 'firstGene')

########################################################################
######## Adding orthologous genes and orthologous interactions #########
########################################################################

full_annotation['isTgOrthologous'] = full_annotation['TGs'].apply(setting_orthologous_tgs, orto_gene=orthologous_gene)
full_annotation = pd.merge(full_annotation, network_info[[0,1]].drop_duplicates(subset=None, keep='first', inplace=False), how = 'left',  left_on='TFname', right_on= 1)
full_annotation = full_annotation.rename(columns={0: 'TFregulonDBid'})
full_annotation = full_annotation.drop(1, axis=1)
all_orthologous_hits = []  # Here we saved all the orthologous interactions
full_annotation['isOrthologousInteraction'] = full_annotation.apply(process_and_log, axis=1)
full_annotation = full_annotation.drop('TFregulonDBid', axis=1)
orthologous_interactions_df = pd.DataFrame(all_orthologous_hits)

full_annotation.to_csv('/space24/PGC/emhernan/3_TRN/output/oTF-TGs-TRN-bgM1.tsv', index=False, sep = '\t')

########################################################################
################## Computing metrics for the paper #####################
########################################################################

TF_TGs = full_annotation.groupby('TFname').apply(summarize_tf).reset_index(drop=True)
TF_TGs.to_csv('/space24/PGC/emhernan/3_TRN/tables/MainTab1.tmp', index=False, sep = '\t')

# Number of TF-TG interactions
tf_tg_num = sum(TF_TGs['number_of_tgs'])
print(f'There are {tf_tg_num} TF-TG interactions')

# Number of average TGs regulated by each TF
tgs_per_tf_aver_num = TF_TGs['number_of_tgs'].mean()
print(f'In avergare each TF is regulating {tgs_per_tf_aver_num} TGs')

# Number of unique TGs
all_tgs = ','.join(full_annotation['TGs'].dropna())
all_tgs_list = [tag.strip() for tag in all_tgs.split(',')]
unique_tgs = set(all_tgs_list)
print(f"There are {len(unique_tgs)} unique TGs")

# Number of orthologous TGs
tgs_orthologous_uniq = setting_orthologous_unqiue_tgs(",".join(unique_tgs), orto_gene=orthologous_gene)
print(f"There are {tgs_orthologous_uniq} orthologous unique TGs")

# Number of orthologous TF-TG interactions
print(f"There are {len(orthologous_interactions_df)} orthologous TF-TG interactions")

# Number of orthologous TF-TG interactions self-regulations and their TF involved
tf_tg_orto = pd.merge(orthologous_interactions_df, trn_first_gene[['TFname', 'TFoldRZlocusTag']].drop_duplicates(subset=None, keep='first', inplace=False), how = 'left', on = 'TFname')
tf_tg_orto['isOrtoSelf'] = np.where(tf_tg_orto['TG_locusTag'] == tf_tg_orto['TFoldRZlocusTag'], 'T', 'F')
print(f"There are {sum(tf_tg_orto['isOrtoSelf'] == 'T')} orthologous TF-TG self-interactions")
TFs = ",".join(tf_tg_orto[tf_tg_orto['isOrtoSelf'] == 'T']['TFname'].unique())
print(f'The orthologous TF-TG self-interactions involve the TFs: {TFs}')