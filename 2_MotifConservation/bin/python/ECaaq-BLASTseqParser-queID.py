import pandas as pd
import os
import re
from collections import defaultdict

path = '/home/emhernan/2_MotifConservation/formated/ECaaqID_temp.tsv'
outpath = '/home/emhernan/2_MotifConservation/formated/ECaaq-queID.tsv'
queID= defaultdict(int)
ident = 1


with open(path, 'r') as file:
    for line in file:
        if(re.match("gnl", line)):
            line = line.replace('\n','')
            key = str(line.split('\t')[0])
        if(re.match("Score", line)):
            queID[key]+=1

with open(outpath, 'w') as ofile:
    for key in queID:
        s = (key + '\n') * queID[key]
        ofile.write(s)