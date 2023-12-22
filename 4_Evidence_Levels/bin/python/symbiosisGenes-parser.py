import pandas as pd
import os
import re
from collections import defaultdict

path = '/home/emhernan/4_Evidence_Levels/annotation/symbiosisGenes.tmp'
dic={'geneName':[], 'locusTag':[], 'proteinID':[]}
newGen = False


with open (path,'r') as f:
    for line in f:
        #print(line)
        if(re. findall("CDS=.*", line)):
            newGen = True
        if(newGen):
            if(re. findall("/gene=", line)):
                pattern = re.search(r"\/gene=\"(.*)\"", line)
                dic['geneName'].append(pattern.group(1))
            if(re. findall("/locus_tag=", line)):
                pattern = re.search(r"\/locus_tag=\"(.*)\"", line)
                dic['locusTag'].append(pattern.group(1))
            if(re. findall("/protein_id=", line)):
                pattern = re.search(r"\/protein_id=\"(.*)\"", line)
                dic['proteinID'].append(pattern.group(1))
                newGen = False



df = pd.DataFrame(data = dic)

df.to_csv('/home/emhernan/4_Evidence_Levels/tables/SuppTable2.tsv', sep="\t", index = False)