import pandas as pd
import os
import re
from collections import defaultdict

path = '/home/emhernan/1_BBH_TFs/motifsInfo/motifs_seq_relation.tsv'
dic={'geneName':[], 'locusTag':[], 'motifDesc':[], 'start':[], 'end':[]}


with open (path,'r') as f:
    for line in f:
        geneName = line.split('\t')[0]
        locusTag = line.split('\t')[1]
        motifDesc=line.split('\t')[2]
        
        if(len(re. findall("\/\/", motifDesc)) == 0):
            
            pattern = re.search(r"(\w+\-*\s*\w+\-*\w+\-*\w+)\s+(\d*)\-*\>*(\d*)", motifDesc)
            motifDescD = pattern.group(1)
            start = pattern.group(2)
            end = pattern.group(3)
            
            dic['geneName'].append(geneName)
            dic['locusTag'].append(locusTag)
            dic['motifDesc'].append(motifDescD)
            dic['start'].append(start)
            dic['end'].append(end)
            
        else:
            for i in range(len(re. findall("\/\/", motifDesc))+1):
                motifTemp = motifDesc.split("//")[i]
                pattern = re.search(r"\s*(\w+\-*\s*\w*\s*\-*\w*\s*\-*\w*\s*\w*\s*\w*\;*\-*\w*\'*)\s+(\d*)\-*\>*(\d*)", motifTemp)
                motifDescD = pattern.group(1)
                start = pattern.group(2)
                end = pattern.group(3)
                dic['geneName'].append(geneName)
                dic['locusTag'].append(locusTag)
                dic['motifDesc'].append(motifDescD)
                dic['start'].append(start)
                dic['end'].append(end)


df = pd.DataFrame(data = dic)
df.to_csv('/home/emhernan/1_BBH_TFs/motifsInfo/motifs_seq_relation_v2.tsv', sep="\t", index = False)