import pandas as pd
import os
import re

path = '/home/emhernan/2_MotifConservation/formated/ECaaq_RZaadb_blastP_seq_formated.tsv'
outpath = '/home/emhernan/2_MotifConservation/formated/ECaaq_RZaadb_blastP_seq_formated_v2.tsv'
cont = 0

with open (path,'r') as f:
    with open(outpath, 'w') as outf:
        for line in f:
            line = line.replace('\n','')
            if(line.split('\t')[3]):
                start = str(line.split('\t')[0])
                end = str(line.split('\t')[1])
                seq = str(line.split('\t')[2])
                seqIdent = str(line.split('\t')[3])
                indexes = [match.start() for match in re.finditer('-', seq)]
                if(indexes):                 
                    for i in indexes:
                        seq = seq.replace('-','', i)
                        seqIdent = list(seqIdent)
                        seqIdent[i-cont] = ''
                        seqIdent = "".join(seqIdent)
                        cont+=1
                    outf.write(start+'\t'+end+'\t'+seq+'\t'+seqIdent+'\n')
                    cont = 0
                else:
                    outf.write(line+'\n')