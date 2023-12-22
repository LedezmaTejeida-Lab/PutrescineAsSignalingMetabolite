indir=/home/emhernan/4_Evidence_Levels/
bin=/home/emhernan/4_Evidence_Levels/bin/
annotation=/home/emhernan/4_Evidence_Levels/annotation/
tables=/home/emhernan/4_Evidence_Levels/tables/
png=/home/emhernan/4_Evidence_Levels/png/
input=/home/emhernan/4_Evidence_Levels/input/

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
