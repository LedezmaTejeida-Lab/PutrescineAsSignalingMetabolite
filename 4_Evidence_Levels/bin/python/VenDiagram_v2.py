import matplotlib as plt
from matplotlib_venn import venn3, venn3_circles
import numpy as np




GeneExpresionPath = "/home/emhernan/4_Evidence_Levels/input/GeneExpression.tsv"
RegulationPath = "/home/emhernan/4_Evidence_Levels/input/oTF-TG-TRN-10-5-bgM1.tsv"
FunctionPath = "/home/emhernan/4_Evidence_Levels/tables/SuppTable2.tsv"


GeneExpresion = []
TranscriptionalRegulation = []
Function = []


with open (GeneExpresionPath,'r') as f:
    for line in f:
        line = line.replace('\n','')
        LT = line.split('\t')[2]
        GeneExpresion.append(LT)

        
with open (RegulationPath,'r') as f:
    for line in f:
        line = line.replace('\n','')
        LT = line.split('\t')[12]
        TranscriptionalRegulation.append(LT)

with open (FunctionPath,'r') as f:
    for line in f:
        line = line.replace('\n','')
        LT = line.split('\t')[1]
        Function.append(LT)

GeneExpresion.remove('LocusTag') 
TranscriptionalRegulation.remove('TG_oldRZlocusTag')
Function.remove('locusTag')

GeneExpresion = set(GeneExpresion)
TranscriptionalRegulation = set(TranscriptionalRegulation)
Function = set(Function)
GeneExpresion.remove('')


fig = plt.pyplot.figure(figsize=(4, 6), dpi=300, constrained_layout=True)
venn_diagram = venn3(subsets = [GeneExpresion, TranscriptionalRegulation, Function], 
      set_labels =('Gene\nExpresion','Transcriptional\nRegulation', 'Function'))



for t in venn_diagram.set_labels: t.set_fontsize(10)
for t in venn_diagram.subset_labels: t.set_fontsize(8)
venn_diagram.get_label_by_id("100").set_x(-0.25)
venn_diagram.get_label_by_id("111").set_y(-0.25)
venn_diagram.get_label_by_id("101").set_y(-0.30)
venn_diagram.get_label_by_id("011").set_y(-0.27)
venn_diagram.get_label_by_id("011").set_x(0.21)
venn_diagram.get_label_by_id("001").set_y(-0.35)
venn_diagram.get_label_by_id("010").set_x(0.45)


plt.pyplot.text(0.5,-0.5, '${n=1661}$',fontsize = 10)
fig.patch.set_alpha(0)
plt.pyplot.savefig('/home/emhernan/4_Evidence_Levels/png/MainF5.png', dpi = 300, bbox_inches='tight')
plt.pyplot.savefig('/home/emhernan/4_Evidence_Levels/png/MainF5.pdf',  format='pdf', dpi = 300, bbox_inches='tight')
plt.pyplot.show()