import pandas as pd
import os
import re
from collections import defaultdict

path = '/home/emhernan/2_MotifConservation/BLASTresultsSeq/ECaaq_RZaadb_blastP_b1_m0_seq.txt'
outpath = '/home/emhernan/2_MotifConservation/formated/ECaaq-subjID.tsv'
subjID= defaultdict(int)
ident = 1

with open(path, 'r') as file:
    for line in file:
        if(re.match(">", line)):
            line = line.replace('\n','')
            line = line.replace('>','')
            line = line.replace(' ','\t')
            key = str(line.split('\t')[0]) + '_' + str(ident)
            ident+=1
        if(re.match("\s+Score\s*= \s*", line)):
            subjID[key]+=1

with open(outpath, 'w') as ofile:
    for key in subjID:
        s = (key.split('_')[0] + '\n') * subjID[key]
        ofile.write(s)