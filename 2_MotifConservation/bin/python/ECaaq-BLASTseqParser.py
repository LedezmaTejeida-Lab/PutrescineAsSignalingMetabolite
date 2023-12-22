import pandas as pd
import os
import re
from collections import defaultdict


# Defining variables

path = '/home/emhernan/2_MotifConservation/formated/ECaaq_RZaadb_blastP_seq.tab'
cont = 0
dic={'start':[], 'end':[], 'seq':[], 'seqIdent':[]}
tmpS = -1
tmpE = -1
tmpSIdent = -1
tmpEIdent = -1
tmpseq = ''
tmpseqIdent = ''
NewSeq = 0

with open (path,'r') as f:
    for line in f:
        if(line.split('\t')[0] == "Score"):
            NewSeq = 1
        if(line.split('\t')[0] and line.split('\t')[0] != "Score"):
            line = line.replace('\n','')
            start = int(line.split('\t')[0])
            end = int(line.split('\t')[2])
            seq = line.split('\t')[1]
            if tmpE+1 == start and NewSeq != 1:
                tmpseq = tmpseq + seq
                tmpE = end
                
            else:
                dic['start'].append(tmpS)
                dic['end'].append(tmpE)
                dic['seq'].append(tmpseq)
                tmpS = start
                tmpE = end
                tmpseq = seq
            
        if(line.split('\t')[0] == ''):
            line = line.replace('\n','')
            seqIdent = str(line.split('\t')[1])                
            dif = int(len(seq)) - int(len(seqIdent))
            seqIdent = ' '* dif + str(seqIdent)  
            if tmpEIdent+1 == start and NewSeq != 1:
                tmpseqIdent = tmpseqIdent + seqIdent
                tmpEIdent = end
                
            else:
                dic['seqIdent'].append(tmpseqIdent)
                tmpseqIdent = seqIdent
                tmpSIdent = start
                tmpEIdent = end
                NewSeq = 0
    dic['start'].append(tmpS)
    dic['end'].append(tmpE)
    dic['seq'].append(tmpseq)
    dic['seqIdent'].append(tmpseqIdent)

df = pd.DataFrame(data = dic)

# saving as tsv file
df.to_csv('/home/emhernan/2_MotifConservation/formated/ECaaq_RZaadb_blastP_seq_formated.tsv', sep="\t", index = False)