import pandas as pd
import os
import re
from collections import defaultdict


# Defining variables
path = '/home/emhernan/2_MotifConservation/formated/RZaaq_ECaadb_blastP_seq.tab'
cont = 0
dic={'start':[], 'end':[], 'seq':[], 'seqIdent':[]}
tmpE = -1
tmpseq = ''
NewSeq = 0
tmpseqIdent = ''

with open (path,'r') as f:
    for line in f:
        if(line.split('\t')[0] == "Score"):
            NewSeq = 1
            dic['seqIdent'].append(tmpseqIdent)
            dic['end'].append(tmpE)
            dic['seq'].append(tmpseq)
            tmpseqIdent = ''
            #print(line)
        if(line.split('\t')[0] == '' and line.split('\t')[1]):
            line = line.replace('\n','')
            # print(line)
            seqIdent = str(line.split('\t')[1])               
        if(line.split('\t')[0] and line.split('\t')[0] != "Score"):
            line = line.replace('\n','')
            # print(line)
            start = int(line.split('\t')[0])
            end = int(line.split('\t')[2])
            seq = line.split('\t')[1]
            dif = int(len(seq)) - int(len(seqIdent))
            seqIdent = ' '* dif + str(seqIdent)
            tmpseqIdent = tmpseqIdent + seqIdent  
            if tmpE+1 == start and NewSeq == 0:
                tmpseq = tmpseq + seq
                tmpE = end
                
            else:
                dic['start'].append(start)
                tmpE = end
                tmpseq = seq
                NewSeq = 0
            

    dic['end'].append(tmpE)
    dic['seq'].append(tmpseq)
    dic['seqIdent'].append(tmpseqIdent)

df = pd.DataFrame(data = dic)

# saving as tsv file
df.to_csv('/home/emhernan/2_MotifConservation/formated/RZaaq_ECaadb_blastP_seq_formated.tsv', sep="\t", index = False)
