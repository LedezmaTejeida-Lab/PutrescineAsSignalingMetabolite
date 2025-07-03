import pandas as pd
import os
import re
from collections import defaultdict

path = '/space24/PGC/emhernan/3_TRN/output/temporalRelationship.tmp'
dic={'TG_RZlocusTag':[], 'TG_oldRZlocusTag':[]}



with open (path,'r') as f:
    for line in f:
        if(re.findall("/locus_tag=", line)):
            pattern = re.search(r"\/locus_tag=\"(.*)\"", line)
            dic['TG_RZlocusTag'].append(pattern.group(1))
        if(re.findall("/old_locus_tag=", line)):
            pattern = re.search(r"\/old_locus_tag=\"(.*)\"", line)
            dic['TG_oldRZlocusTag'].append(pattern.group(1))
        if (len(re.findall("/locus_tag=", line)) == 0 and len(re.findall("/old_locus_tag=", line)) == 0):
            dic['TG_oldRZlocusTag'].append("NA")

df = pd.DataFrame(data = dic)
df2 = df.drop_duplicates()
df2.to_csv('/space24/PGC/emhernan/3_TRN/output/p4_tmp', sep="\t", index = False)
