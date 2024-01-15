# Name:
#  4_Evidence_Levels.sh
# Author:
#  Hernandez-Benitez Ericka Montserrat
# Version
#  v1.5
# 
# Description
#	The script performs the following analyses:
#		1) Identify nod, fix and nif genes from the annotation file.  
#		2) Plots: Venn Diagram (Main Fig 5), Heatmap (Supp Fig 4) and Supp Table 2
# sh 4_Evidence_Levels.sh


# Additional notes:
# 	The script contemplates the following structure in directories and thus
#	the script will NOT generate them.

#.
#├── annotation
#│   └── AHJU02.1.gbff
#├── bin
#│   ├── bash
#│   │   └── 4_Evidence_Levels.sh
#│   ├── python
#│   │   ├── symbiosisGenes-parser.py
#│   │   └── VenDiagram_v2.py
#│   └── R
#│       └── NodulationHeatmap.R
#└── input
#   ├── GeneExpression.tsv
#   └── oTF-TG-TRN-10-5-bgM1.tsv


indir=/home/emhernan/4_Evidence_Levels/
bin=/home/emhernan/4_Evidence_Levels/bin/ # Must be created with all the binary files
annotation=/home/emhernan/4_Evidence_Levels/annotation/ # Must be created with all the annotation files
tables=/home/emhernan/4_Evidence_Levels/tables/
png=/home/emhernan/4_Evidence_Levels/png/
input=/home/emhernan/4_Evidence_Levels/input/ # Must be created

cd $indir


# 1) Generating Supplementary Table 2

grep -A8 -wE "/gene=\"nod[A-Z]*[0-9]*[a-z]*\"|/gene=\"nif[A-Z]*[0-9]*[a-z]*\"|/gene=\"fix[A-Z]*[0-9]*[a-z]*\"|\/gene=\"nolO\"|\/gene=\"nfeD\"|\/gene=\"fnrN[a-z]*\"|\/gene=\"nolL\"" $annotation'AHJU02.1.gbff'  | sed 's/                     //g'| sed 's/     //g' | sed 's/   /=/g' | sed 's/\-\-//g' | grep -v "CDS=101404..101847" > $annotation'symbiosisGenes.tmp'

mkdir $tables
python $bin'python/symbiosisGenes-parser.py' && rm $annotation'symbiosisGenes.tmp'


# 2) Generating the intersection plot

mkdir $png
cd $input
python $bin'python/VenDiagram_v2.py'


# 3) Generating Supplementary Figure 4

Rscript $bin'R/NodulationHeatmap.R'
