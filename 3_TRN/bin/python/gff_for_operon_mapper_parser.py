# Loading libraries
import re
import pandas as pd
import numpy as np

# Loading gff file
gff_file = "/space24/PGC/emhernan/3_TRN/operon_mapper/input/GCF_000268285.2_RPHCH2410v2_genomic.gff"
gff = pd.read_csv(gff_file, sep="\t", comment = "#", header = None)
gff = gff[~((gff[1] == 'RefSeq') & (gff[2] == 'region'))]

## Setting all my functions
def extract_geneParentLocusTag(info):
    pattern =  r"ID=(gene\d+)\;Name.*locus_tag=(RPHASCH2410\_.*\d+)\;old_locus_tag=(RPHASCH2410\_.*\d+)"
    pattern_pseudo = r"ID=(gene\d+)\;Name.*locus_tag=(RPHASCH2410\_RS\d+)"
    
    match = re.search(pattern, info)
    if match:
        return pd.Series([match.group(1), match.group(2), match.group(3)])
    
    match = re.search(pattern_pseudo, info)
    
    if match:
        return pd.Series([match.group(1), match.group(2), ""])
    
    return pd.Series(["", "", ""])  # fallback if not found pattern

def extract_cds_geneParent(info):
    pattern = r"ID=(cds\d+);Parent=(gene\d+);Dbxref=.*?protein_id=(WP_\d+\.1);"
    pattern_pseudo = r"ID=(cds\d+);Parent=(gene\d+);.*?pseudo=true;"
    
    match = re.search(pattern, info)
    if match:
        return pd.Series([match.group(1), match.group(2), match.group(3)])
    
    match = re.search(pattern_pseudo, info)
    
    if match:
        return pd.Series([match.group(1), match.group(2), ""])
    
    return pd.Series(["", "", ""])  # fallback if not found pattern


def add_ID_column(df):
    df['ID'] = np.where(
        df['old_locus_tag'] != '',
        'ID=' + df['old_locus_tag'] + ';',
        'ID=' + df['locus_tag'] + ';'
    )
    return df

def extract_rna_geneParent(info):
    pattern = r"ID=(rna\d+)\;Parent=(gene\d+)\;.*"
    match = re.search(pattern, info)
    if match:
        return pd.Series([match.group(1), match.group(2)])
    return pd.Series(["", ""])  # fallback if not found pattern


# From ParentID to locus_tag
annotation_parentGene_locusTag = gff[(gff[1] == "RefSeq") & (gff[2] != "region")].copy()
annotation_parentGene_locusTag = annotation_parentGene_locusTag[[3,4,8]]
annotation_parentGene_locusTag[['geneParent', 'locus_tag', 'old_locus_tag']] = annotation_parentGene_locusTag[8].apply(extract_geneParentLocusTag)
annotation_parentGene_locusTag["old_locus_tag"] = annotation_parentGene_locusTag["old_locus_tag"].astype(str).str.replace(";partial.*", "", regex=True)

# From cds to ParentID to Locustag
annotation_cds_parentID = gff[(gff[1] == "Protein Homology") | (gff[1] == "GeneMarkS+")].copy() # 6539 rows 9 columns
annotation_cds_parentID[['cds', 'geneParent', 'proteinID']] = annotation_cds_parentID[8].apply(extract_cds_geneParent)
cds_locus_Tag = pd.merge(annotation_cds_parentID, annotation_parentGene_locusTag, how = 'left', on = ['geneParent',3,4])
cds_locus_Tag = add_ID_column(cds_locus_Tag)
cds_locus_Tag = cds_locus_Tag[[0,1,2,3,4,5,6,7,'ID']]

# From tRNA to ParentID to Locustag
annotation_trna_parentID = gff[gff[2] == "tRNA"].copy() # 51
annotation_trna_parentID[['trna', 'geneParent']] = annotation_trna_parentID[8].apply(extract_rna_geneParent)
trna_locus_Tag = pd.merge(annotation_trna_parentID, annotation_parentGene_locusTag, how = 'left', on = ['geneParent',3,4])
trna_locus_Tag = add_ID_column(trna_locus_Tag)
trna_locus_Tag = trna_locus_Tag[[0,1,2,3,4,5,6,7,'ID']]

# From tRNA to ParentID to Locustag
annotation_rrna_parentID = gff[gff[2] == "rRNA"].copy() # 51
annotation_rrna_parentID[['rrna', 'geneParent']] = annotation_rrna_parentID[8].apply(extract_rna_geneParent)
rrna_locus_Tag = pd.merge(annotation_rrna_parentID, annotation_parentGene_locusTag, how = 'left', on = ['geneParent',3,4])
rrna_locus_Tag = add_ID_column(rrna_locus_Tag)
rrna_locus_Tag = rrna_locus_Tag[[0,1,2,3,4,5,6,7,'ID']]

gff_formatted = pd.concat([cds_locus_Tag, trna_locus_Tag, rrna_locus_Tag], ignore_index=True)
gff_formatted.to_csv('/space24/PGC/emhernan/3_TRN/operon_mapper/input/gff_formatted.tsv', index=False, sep = "\t", header = False)