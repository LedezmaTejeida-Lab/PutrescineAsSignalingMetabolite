# Loading libraries
import re
import pandas as pd

## Setting file paths
operons_file = "/space24/PGC/emhernan/3_TRN/operon_mapper/output/list_of_operons_546291"
operons = pd.read_csv(operons_file, sep="\t")

## Extracting metrics
operons["Operon"] = operons["Operon"].ffill().astype(int)
operons = operons.loc[~operons['IdGene'].isna()]
operon_genes = operons['Operon'].value_counts()

operon_genes = operons['Operon'].value_counts()
print("Average number of genes per operon " + str(operon_genes.mean()))

## Estracting heads of operons
operons['PosLeft'] = pd.to_numeric(operons['PosLeft'], errors='coerce')

# Split by Strand
strand_plus = operons[operons['Strand'] == '+']
strand_minus = operons[operons['Strand'] == '-']

# Obtaining the head of operons by strans
head_plus = strand_plus.loc[strand_plus.groupby('Operon')['PosLeft'].idxmin()]
head_minus = strand_minus.loc[strand_minus.groupby('Operon')['PosLeft'].idxmax()]

# Binding heads of operons
operon_heads = pd.concat([head_plus, head_minus]).sort_values('Operon')

# Retrieving only OperonId and GeneId
operon_heads = operon_heads[['Operon', 'IdGene']]

# Saving in a tsv file
operon_heads.to_csv('/space24/PGC/emhernan/3_TRN/operon_mapper/output/head_operons_locus_tag_gff_formatted.tsv', index=False, sep = "\t", header = False)