import pandas as pd
import os
import re
from collections import defaultdict

path = '/home/emhernan/1_BBH_TFs/tables/TF_Orthologous_RZannotation_table.tsv'
dic={'product':[], 'protein_id':[]}


with open (path,'r') as f:
    for line in f:
        line =  line.lstrip()
        if(re.findall("/product=", line)):
            pattern = re.search(r"\/product=\"(.*)", line)
            annotation = pattern.group(1)
            annotation = annotation.replace("\"","")
        if not re.match(r"^/", line) and not re.match(r"--", line):
            pattern2 = re.search(r"(.*)", line)
            annotation = annotation + " " + pattern2.group(1)
            
        if(re.findall("protein_id=", line)):
            annotation = annotation.replace("\"","")
            dic['product'].append(annotation)
            pattern3 = re.search(r"\/protein_id=\"(.*)\"", line)
            p_id = pattern3.group(1)
            dic['protein_id'].append(p_id)

df = pd.DataFrame(data = dic)
df.to_csv('/home/emhernan/1_BBH_TFs/tables/TF_Orthologous_RZannotation_table_v2.tsv', sep="\t", index = False)