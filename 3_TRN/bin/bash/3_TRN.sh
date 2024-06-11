# Name:
#  3_TRN.sh
# Author:
#  Hernandez-Benitez Ericka Montserrat
# Version
#  v1.3
# 
# Description
#	The script performs the following analyses:
#		1) Formating PSSM matrixes according to RSAT input requirements. 
#		2) Filtering the PSSM of TF with a DNA-binding region motif with identity >= 30 and coverage = 100 in both directions.
#		3) Retrieve upstream sequences from all R. phaseoli genes and compute matrix scan with the filtered PSSMs. 
#		4) Create the file oTF-TG-TRN-10-5-bgM1.tsv
#		5) Plots: Bar plot TF vs. TG ortholous and non-orthologous (Main Fig 4) and Main Table 1.
# sh 3_TRN.sh


# Additional notes:
# 	The script contemplates the following structure in directories and thus
#	the script will NOT generate them.

#.
#├── annotation
#│   ├── GCA_000268285.2_genomic.gbff
#│   └── GCF_000268285.2_genomic.gbff
#├── bin
#│   ├── bash
#│   │   └── 3_TRN.sh
#│   ├── python
#│   │   └── TG-LocusTag-parsey.py
#│   └── R
#│       ├── DNAbindingFilterTFs100Cov30Ident.R #change
#│       ├── merge_MatrixScantmpFiles.R
#│       ├── oTF-TG-NetworkGenerator.R
#│       ├── section3-MainFigure4.R
#│       └── section3-MainTable1.R
#├── matrixes
#│   ├── PSSM-Dataset-v4.0-ed.txt
#│   ├── PSSM-Dataset-v4.0.txt
#│   ├── RhizobiumphaseoliOrthologousCutOffIdent30Cov100Matrices.tab
#│   └── RhizobiumphaseoliOrthologousMatrices.tab
#├── networkInfo
#│   └── network_tf_gene.txt
#└── OrthologousTFInfo
#   ├── ECaaq_oTFs_DNAbinding_motifs_info.tsv
#   ├── ECaaq_oTFs_motifs_info.tsv
#   ├── Orthologous_table.tsv
#   ├── RZaaq_oTFs_DNAbinding_motifs_info.tsv
#   ├── RZaaq_oTFs_motifs_info.tsv
#   └── TF_Orthologous_table.tsv

#########################################################################################################
########################################## Working directory ###########################################
#########################################################################################################


cd /space24/PGC/emhernan/3_TRN/


#########################################################################################################
########################################### Saving directories ##########################################
#########################################################################################################
indir=/space24/PGC/emhernan/3_TRN/
bin=/space24/PGC/emhernan/3_TRN/bin/ # Must be created with all the binaries files
matrixes=/space24/PGC/emhernan/3_TRN/matrixes/ # Must be created with Regulondb matrixes (PSSMs)
networkInfo=/space24/PGC/emhernan/3_TRN/networkInfo/ # Must be created with network information
OrthologousTFInfo=/space24/PGC/emhernan/3_TRN/OrthologousTFInfo/ # Must be created with all orthologous information
annotation=/space24/PGC/emhernan/3_TRN/annotation/ # Must be created with all annotation information
sequences=/space24/PGC/emhernan/3_TRN/sequences/
output=/space24/PGC/emhernan/3_TRN/output/
png=/space24/PGC/emhernan/3_TRN/png/
tables=/space24/PGC/emhernan/3_TRN/tables/
bgfilePath=/space23/rsat/rsat/public_html/data/genomes/Rhizobium_phaseoli_GCF_000268285.2_RPHCH2410v2/oligo-frequencies/2nt_upstream-noorf_Rhizobium_phaseoli_GCF_000268285.2_RPHCH2410v2-ovlp-1str.freq.gz


cd $indir

# 1)  Matrixes file iles directory
grep -v "#" $matrixes'PSSM-Dataset-v4.0.txt' | perl -nae 'if(/Transcription Factor Name: (\w+)/) {print ";$1\tTranscription Factor\n"} else {if(/([ACG])(\t[\d+\t+]*)/) {print "$1$2\n"} else {if(/(T)(\t[\d+\t+]*)/) {print "$1$2\n//\n"}}}' > $matrixes'PSSM-Dataset-v4.0-ed.txt'

# 2) Filtering TFs with identity >= 40 and coverage = 100 in their DNA-binding motif
Rscript --vanilla $bin'R/DNAbindingFilterTFs100Cov40Ident.R' $OrthologousTFInfo 'TF_Orthologous_table_filter100Cov40Iden.tsv'

# 3) Filtering matrixes to 17 TFs with PSSM and with identity >= 40 and coverage = 100 in their DNA-binding motif
grep -A 4 -iwf <(grep -iwf <(grep "Transcription Factor Name:" $matrixes'PSSM-Dataset-v4.0.txt'  | cut -f2 -d ":" | sed 's/ //g' | sort | uniq) <(cat $OrthologousTFInfo'TF_Orthologous_table_filter100Cov40Iden.tsv') | cut -f2 | sed 's/glnG/ntrC/' | sed 's/ttdR/dan/g') <(grep -v "#" $matrixes'PSSM-Dataset-v4.0.txt' | perl -nae 'if(/Transcription Factor Name: (\w+)/) {print "; $1\tTranscription Factor\n"} else {if(/([ACG])(\t[\d+\t+]*)/) {print "$1 |$2\n"} else {if(/(T)(\t[\d+\t+]*)/) {print "$1 |$2\n\/\/"}}}')  | grep -v "\-\-" |sed  's/\/\/;/\/\/\n;/g' > $matrixes'RhizobiumphaseoliOrthologousCutOffIdent40Cov100Matrices.tab'


#########################################################################################################
##################################### Peforming matrix-scan analysis ####################################
#########################################################################################################

# 1) Retrive sequences
# To work with a genome for which no mRNA is available, choose ‘CDS’ as feature type and upstream sequences will be retrieved relative to the start codon
# Change organism for Rhizobium_phaseoli_GCF_000268285.2_RPHCH2410v2
mkdir $sequences
rsat retrieve-seq -org Rhizobium_phaseoli_GCF_000268285.2_RPHCH2410v2  -feattype CDS -type upstream -format fasta -label id,name -from -400 -to +50 -noorf -all > $sequences'retrieve-seq-Rhizobium_phaseoli_Ch24-10-50.fna'

# 2) Convert file to fasta format with Convert-seq
rsat convert-seq  -i $sequences'retrieve-seq-Rhizobium_phaseoli_Ch24-10-50.fna' -from  fasta -to fasta  -mask non-dna -o $sequences'retrieve-seq-Rhizobium_phaseoli_Ch24-10-50.fnaed'

# 3) Matrix-scan with 17 TF PSSMs
mkdir $output
rsat matrix-scan -v 1 -matrix_format tab -m $matrixes'RhizobiumphaseoliOrthologousCutOffIdent40Cov100Matrices.tab' -consensus_name  -pseudo 1 -decimals 1 -2str -origin end -bgfile $bgfilePath -bg_pseudo 0.01 -return limits -return sites -return pval -return rank -lth score 1  -uth pval 1e-5  -i $sequences'retrieve-seq-Rhizobium_phaseoli_Ch24-10-50.fnaed' -seq_format fasta -n score > $output'scan-result-10-5-50-bgM1.txt' &

# 4) Feature map
# Note: -A 18  = 17 matrixes +1 
paste <(grep -A 18 -wE "Matches" $output'scan-result-10-5-50-bgM1.txt' | cut -f3 | grep -v ";" | grep -v "name" |sed 's/ //g') <(grep ";" $matrixes'RhizobiumphaseoliOrthologousCutOffIdent40Cov100Matrices.tab' | cut -f2 -d " " | cut -f1)  > $output'relation_tmp'
grep -wf <(grep -v ";" $output'scan-result-10-5-50-bgM1.txt' | grep -v "#" | cut -f1-6 | grep -w "site" | cut -f1 | sort | uniq -c | sort -nr | cut -f2 -d "|" | sed -n '1,15p')  $output'scan-result-10-5-50-bgM1.txt' | grep -v ";" | cut -f1-6 | grep -w "site" > $output'scan-result-10-5-50-bgM1.tmp'
Rscript --vanilla $bin'R/merge_MatrixScantmpFiles.R' $output'scan-result-10-5-50-bgM1.tmp' $output'relation_tmp' $output'scan-result-10-5-50-bgM1.tmp2'

rsat feature-map  -htmap  -scorethick  -legend  -scalebar  -horizontal  -title 'Motifs discovered in upstream sequences in Rhizobium phaseoli Individual sites (10-5)' -from auto -to auto -origin 0 -mlen 500 -mapthick 25 -mspacing 2 -format png -bgcolor 220,220,220 -i $output'scan-result-10-5-50-bgM1.tmp2'  -o $png'scan-result-10-5-50-bgM1.png' > $output'scan-result-10-5-50-bgM1.html'

cd $output && rm *tmp*


# Important notes: The first 6 columns must be cut off from the scan-matrix output.
##################################################################################
##################################################################################

# Feature map format
# In the standard input format, each feature is represented by a single line, which must contain the following information:

#    - map name (eg: gene name),
#    - feature type (site, ORF),
#    - identifier(ex: GATA_box, Abf1_site)
#    - strand (D for Direct, R for Reverse),
#    - start position (may be negative)
#    - end position (may be negative)
#    - (optional) sequence (ex: AAGATAAGCG)
#    - (optional) score (a real value)

# Fields must be provided in this order, separated by tabs.
# Feature map also imports result files from some other programs or databases: 
##################################################################################
##################################################################################


# The meaning of the numbers
# 10-5-50

#    - 10-5 : treshold 10^-5
#    - 50 : The retrieving sequence process was made from -400 to +50 
## We only considered 10^-5  treshold without CRE since it is unuseful. 
## From our 75 oTF only 35 have an experimental PSSM
##################################################################################
##################################################################################



#########################################################################################################
########################################### Assembling the TRN ##########################################
#########################################################################################################

# 1) Creating temporary files to generate the oTF-TGInteractions.tsv file
# Important notes:
# Change ttdR to dan
# Change glnG to ntrC

# 1.1.1) Generating temporary file TG_RZlocusTag, TG_RZname, typeOfMatch, TF_consensusMotif
grep -v ";" $output'scan-result-10-5-50-bgM1.txt' | cut -f1-3 | grep -v "limit" | sed 's/|/\t/' | grep -v "#seq_id" | sort | uniq > $output'p1_tmp'

# 1.1.2) Generating temporary file TF_consensusMotif, TF_name
paste <(grep -A 18 "Matches per matrix" $output'scan-result-10-5-50-bgM1.txt' |  cut -f3 | grep -v ";" | grep -v "name" |sed 's/ //g') <(grep ";" $matrixes'RhizobiumphaseoliOrthologousCutOffIdent40Cov100Matrices.tab'  |cut -f2 -d " " | cut -f1) > $output'p2_tmp'

# 1.1.3) Generating temporary file TF_locusTag, TF_name, RegulondbID, TF_RZlocusTag, TF_RZoldLocusTag
grep -iwf <(cut -f2 $output'p2_tmp') <(cut -f3,4,9,11,12 $OrthologousTFInfo'TF_Orthologous_table.tsv' | sort | uniq | sed 's/ttdR/dan/g'| sed 's/glnG/ntrC/g') > $output'p3_tmp'

# 1.1.4) Generating temporary file with orthologus TG: TG_locusTag, TG_name, TG_regulonDBid, TG_RZlocusTag
# Double check that grep -wf <(cut -f2 $output'p2_tmp') <(grep -v "#" $networkInfo'network_tf_gene.txt' | cut -f1,2,3,4 | sort | uniq) | cut -f2 | sort | uniq | wc == 17; if not check synonyms
grep -wf <(grep -wf <(cut -f2 $output'p2_tmp') <(grep -v "#" $networkInfo'network_tf_gene.txt' | cut -f1,2,3,4 | sort | uniq) | cut -f4 | sort  |uniq) <(cut -f3,4,5,6,7,8,9,11,12 $OrthologousTFInfo'Orthologous_table.tsv' | sort | uniq)  | cut -f1,2,7,9 > $output'p4_tmp'

# 1.1.5) Generating relationship TG_RZlocusTag TG_oldRZlocusTag
grep  -B1 -wf  <(cut -f1 p1_tmp | sort | uniq) $annotation'GCF_000268285.2_genomic.gbff' | sed 's/                     //' | sed 's/--//'| sed '/^$/d' > $output'temporalRelationship.tmp'
python $bin'python/TG-LocusTag-parsey.py' && rm $output'temporalRelationship.tmp'


# 1.1.6) Generating temporary file TF_RegunlobdbID, TG_RegunlobdbID, TG_name, TFname_TGRegunlondbId
# Double check the grep -wf <(cut -f2 $output'p2_tmp') <(grep -v "#" $networkInfo'network_tf_gene.txt' | cut -f1,2,3,4 | sort | uniq) | cut -f2 | sort | uniq | wc == 17; if not check synonyms

grep -wf <(cut -f9 $OrthologousTFInfo'Orthologous_table.tsv' | sort | uniq) <(paste <(grep -wf <(cut -f2 $output'p2_tmp') <(grep -v "#" $networkInfo'network_tf_gene.txt' | cut -f1,2,3,4 | sort | uniq) | cut -f1,3,4) <(grep -wf <(cut -f2 $output'p2_tmp') <(grep -v "#" $networkInfo'network_tf_gene.txt' | cut -f1,2,3,4 | sort | uniq) | cut -f2,3 | sed 's/\t/-/g')) > $output'p5_tmp'
Rscript --vanilla $bin'R/oTF-TG-NetworkGenerator.R' $output'p1_tmp' $output'p2_tmp' $output'p3_tmp' $output'p4_tmp' $output'p5_tmp' $output'p6_tmp' $output'oTF-TG-TRN-10-5-bgM1.tsv' && cd $output && rm *tmp

#########################################################################################################
#################################### Plotting the TRN (Main figure 4) ###################################
#########################################################################################################

mkdir $png
Rscript --vanilla $bin'R/section3-MainFigure4.R'  $output'oTF-TG-TRN-10-5-bgM1.tsv' $OrthologousTFInfo'Orthologous_table.tsv' $png


#########################################################################################################
#################################### Computing table 1 (Main table 1) ###################################
#########################################################################################################


mkdir $tables
Rscript --vanilla $bin'R/section3-MainTable1.R' $output'oTF-TG-TRN-10-5-bgM1.tsv' $OrthologousTFInfo'TF_Orthologous_table.tsv' $OrthologousTFInfo'ECaaq_oTFs_motifs_info.tsv' $OrthologousTFInfo'RZaaq_oTFs_motifs_info.tsv' $tables'MainTab1.tsv'



#########################################################################################################
############################################## Computing control ########################################
#########################################################################################################

# cd 3_TRN/output-tmp

#TFsMotifs=$(grep -wif <(grep -w ";" $matrixes'RhizobiumphaseoliOrthologousCutOffIdent40Cov100Matrices.tab' | cut -f2 -d " " | cut -f1 | sort  |uniq) relation_tmp-2021-11-23 | cut -f1)
#matrixMatch=$(grep -wf <(echo "$TFsMotifs") <(grep -w "site" scan-result-10-5-50-bgM1-35PSSMs-2021-11-23.txt))
#MatrixMotif=$(cat <(grep -A 1 -wE "Matches" scan-result-10-5-50-bgM1-35PSSMs-2021-11-23.txt) <(grep -wf <(echo "$TFsMotifs") <(grep -A 36 -wE "Matches" scan-result-10-5-50-bgM1-35PSSMs-2021-11-23.txt)))
#cat <(echo "$matrixMatch") <(echo "$MatrixMotif") > $output'scan-result-10-5-50-bgM1.txt'
#cat <(echo "$matrixMatch") <(echo "$MatrixMotif") > scan-result-10-5-50-bgM1.txt

#grep -wif <(grep -w ";" $matrixes'RhizobiumphaseoliOrthologousCutOffIdent40Cov100Matrices.tab' | cut -f2 -d " " | cut -f1 | sort  |uniq) relation_tmp-2021-11-23 > $output'relation_tmp'