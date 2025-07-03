# Name:
#  3_TRN.sh
# Author:
#  Hernandez-Benitez Ericka Montserrat
# Version
#  v1.4
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

# .
# ├── annotation
# │   └── GCF_000268285.2_RPHCH2410v2_genomic.gbff
# ├── bin
# │   ├── bash
# │   │   └── 3_TRN.sh
# │   ├── python
# │   │   ├── extending_trn_to_operons.py
# │   │   ├── gff_for_operon_mapper_parser.py
# │   │   ├── operon_mapper_metrics_and_head_of_operons_parser.py
# │   │   └── TG-LocusTag-parsey.py
# │   └── R
# │       ├── DNAbindingFilterTFs100Cov40Ident.R
# │       ├── oTF-TG-NetworkGenerator.R
# │       ├── section3-MainFigure4.R
# │       └── section3-MainTable1.R
# ├── matrixes
# │   ├── PSSM-Dataset-v4.0.txt
# │   └── PSSM-ThreholdsbyTF_v4.0.txt
# ├── networkInfo
# │   └── network_tf_gene.txt
# ├── operon_mapper
# │   ├── input
# │   │   ├── GCF_000268285.2_RPHCH2410v2_genomic.fna
# │   │   └── GCF_000268285.2_RPHCH2410v2_genomic.gff
# │   └── output
# │       ├── functional_descriptions_546291
# │       ├── list_of_operons_546291
# │       ├── operonic_gene_pairs_546291
# │       ├── ORFs_coordinates_546291
# │       ├── predicted_COGs_546291
# │       ├── predicted_orfs_546291
# │       └── predicted_protein_sequences_546291
# └── OrthologousTFInfo
#     ├── ECaaq_oTFs_DNAbinding_motifs_info.tsv
#     ├── ECaaq_oTFs_motifs_info.tsv
#     ├── Orthologous_table.tsv
#     ├── RZaaq_oTFs_DNAbinding_motifs_info.tsv
#     ├── RZaaq_oTFs_motifs_info.tsv
#     └── TF_Orthologous_table.tsv






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
matrixes_split=/space24/PGC/emhernan/3_TRN/matrixes/matrixes_split/ # Must be created with Regulondb matrixes (PSSMs)
networkInfo=/space24/PGC/emhernan/3_TRN/networkInfo/ # Must be created with network information
OrthologousTFInfo=/space24/PGC/emhernan/3_TRN/OrthologousTFInfo/ # Must be created with all orthologous information
annotation=/space24/PGC/emhernan/3_TRN/annotation/ # Must be created with all annotation information
sequences=/space24/PGC/emhernan/3_TRN/sequences/
output=/space24/PGC/emhernan/3_TRN/output/
png=/space24/PGC/emhernan/3_TRN/png/
tables=/space24/PGC/emhernan/3_TRN/tables/
bgfilePath=/space23/rsat/rsat/public_html/data/genomes/Rhizobium_phaseoli_GCF_000268285.2_RPHCH2410v2/oligo-frequencies/2nt_upstream-noorf_Rhizobium_phaseoli_GCF_000268285.2_RPHCH2410v2-ovlp-1str.freq.gz
operonMapper=/space24/PGC/emhernan/3_TRN/operon_mapper/ # Must be created with operon Mapper input files

mkdir -p $tables
mkdir -p $sequences
mkdir -p $matrixes_split
mkdir -p $output
mkdir -p $png

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
rsat retrieve-seq -org Rhizobium_phaseoli_GCF_000268285.2_RPHCH2410v2  -feattype CDS -type upstream -format fasta -label id,name -from -400 -to +50 -noorf -all > $sequences'retrieve-seq-Rhizobium_phaseoli_Ch24-10-50.fna'


# 2) Convert file to fasta format with Convert-seq
rsat convert-seq  -i $sequences'retrieve-seq-Rhizobium_phaseoli_Ch24-10-50.fna' -from  fasta -to fasta  -mask non-dna -o $sequences'retrieve-seq-Rhizobium_phaseoli_Ch24-10-50.fnaed'


# 3) Predict operons to filter the operon leaders from retrive seq output 
# 3.1) Formatting gff file for operon mapper
python3 $bin'python/gff_for_operon_mapper_parser.py'
# 3.2) Provide to operon mapper de gff formated and the genomic.fna

# 4) Filter only the head of operons from retrieve seq output
# 4.1) Extracting head of operond from operon mapper output
python3 $bin'python/operon_mapper_metrics_and_head_of_operons_parser.py'
# Average number of genes per operon 1.6981471950591869
# 4.2) Filter the head of operons from retrieve seq-
awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq=seq $0} END {if (seq) print seq}' $sequences'retrieve-seq-Rhizobium_phaseoli_Ch24-10-50.fnaed' > $sequences'retrieve-seq-Rhizobium_phaseoli_Ch24-10-50.fnaed.tmp'
grep -A1 -wf <(grep -wf <(cut -f2 $operonMapper'/output/head_operons_locus_tag_gff_formatted.tsv') <(cut -f1 $sequences'retrieve-seq-Rhizobium_phaseoli_Ch24-10-50.fnaed.tmp' | sed 's/>//g')) $sequences'retrieve-seq-Rhizobium_phaseoli_Ch24-10-50.fnaed.tmp'| grep -v "\-\-" > $sequences'retrieve-seq-Rhizobium_phaseoli_Ch24-10-50_head_operons.fnaed' && rm $sequences'retrieve-seq-Rhizobium_phaseoli_Ch24-10-50.fnaed.tmp' &
#pause; background process; wc == 3447
# 4.3) Computing the number of pb that scan matrix will scan
grep ">" $sequences'retrieve-seq-Rhizobium_phaseoli_Ch24-10-50_head_operons.fnaed'  | cut -f3 -d ";" | sed 's/ size: //' | awk '{s+=$1} END {print s}'
grep -v ">" $sequences'retrieve-seq-Rhizobium_phaseoli_Ch24-10-50_head_operons.fnaed' | awk '{print length($0)}' | paste -sd+ | bc
# 849,961


# 5)  Running matrix-scan with:
#		a) 17 PSSMs 
#		b) retrieve-seq filtered by head of operons 
#		c) Thresholds established in PSSM-ThreholdsbyTF_v4.0.txt

# 5.1) Splitting PSSM; saving one file per PSSM (17 files are created)
cd $matrixes

for matrix in $(grep ";" "$matrixes"RhizobiumphaseoliOrthologousCutOffIdent40Cov100Matrices.tab  | cut -f2 -d " "| cut -f1)
do
    path=$matrix
    path+="-PSSM-v4.0-ed.txt"
    grep -A5 -w "$matrix" PSSM-Dataset-v4.0-ed.txt > $matrixes_split$path
done


# 5.2) Running matrix scan with the thesholds established in PSSM-ThreholdsbyTF_v4.0.txt

cd $matrixes_split
for matrix in $(ls "$matrixes_split")
do
    TF=$(echo "$matrix" | cut -f1 -d "-")
    Threshold=$(grep -w $TF $matrixes'PSSM-ThreholdsbyTF_v4.0.txt' | cut -f2  | tr -d '[:space:]')
    t="$TF-$Threshold-bgM1.txt"
	rsat matrix-scan -v 1 -matrix_format tab -m "$matrix" -consensus_name -pseudo 1 -decimals 1 -2str -origin end -bgfile $bgfilePath -bg_pseudo 0.01 -return limits -return sites -return pval -uth pval "${Threshold}" -return rank -i $sequences'retrieve-seq-Rhizobium_phaseoli_Ch24-10-50_head_operons.fnaed' -seq_format fasta -n score > $output$t &
    #rsat matrix-scan -v 1 -matrix_format tab -m "$matrix" -consensus_name -pseudo 1 -decimals 1 -2str -origin end -bgfile $bgfilePath -bg_pseudo 0.01 -return limits -return sites -return pval -uth pval "${Threshold}" -return rank -i $sequences'retrieve-seq-Rhizobium_phaseoli_Ch24-10-50_head_operons.fnaed' -seq_format fasta -n score -o $output$t
done

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

############### Checar hasta aquí

#########################################################################################################
########################################### Assembling the TRN ##########################################
#########################################################################################################
# oTF-TGInteractions.tsv file must have the following columns:
# 1) scanType
# 2) TF_RdbID
# 3) TFname
# 4) TF_locusTag
# 5) TF_RZlocusTag
# 6) TF_oldRZlocusTag
# 7) TG_RdbID
# 8) TG_name
# 9) TG_locusTag
# 10) TG_RZlocusTag
# 11) TG_oldRZlocusTag
# 12) TG_RZname
# 13) OrthologousInteration


# 1) Creating temporary files to generate the oTF-TGInteractions.tsv file
# Important notes:
# Change ttdR to dan
# Change glnG to ntrC

# 1.1.1) Generating a temporary file TG_RZoldlocusTag, TG_RZname, typeOfMatch, TF_name


for result in "$output"/*.txt; do
    TF=$(basename "$result" | cut -d'-' -f1)
    
    grep -v ";" "$result" | cut -f1,2 | grep "site" | sed 's/|/\t/' | \
    awk -v tf="$TF" 'BEGIN { OFS="\t" } { print $0, tf }' >> "$output/p1_tmp"
done

sort "$output/p1_tmp" | uniq > "$output/tmp" && mv "$output/tmp" "$output/p1_tmp"
# 1.1.2) Generating temporary file TF_locusTag, TF_name, RegulondbID, TF_RZlocusTag, TF_RZoldLocusTag
grep -iwf <(cut -f4 "$output/p1_tmp") <(cut -f3,4,9,11,12 $OrthologousTFInfo'TF_Orthologous_table.tsv' | sort | uniq | sed 's/ttdR/dan/g'| sed 's/glnG/ntrC/g') > $output'p2_tmp'

# 1.1.3) Generating temporary file with orthologus TG: TG_locusTag, TG_name, TG_regulonDBid, TG_RZoldlocusTag
# Double check that grep -wf <(cut -f4 "$output/p1_tmp") <(grep -v "#" $networkInfo'network_tf_gene.txt' | cut -f1,2,3,4 | sort | uniq) | cut -f2 | sort | uniq | wc == 16; if not check synonyms
grep -wf <(grep -wf <(cut -f4 "$output/p1_tmp") <(grep -v "#" $networkInfo'network_tf_gene.txt' | cut -f1,2,3,4 | sort | uniq) | cut -f4 | sort  |uniq) <(cut -f3,4,5,6,7,8,9,11,12 $OrthologousTFInfo'Orthologous_table.tsv' | sort | uniq)  | cut -f1,2,7,9 > $output'p3_tmp'

# 1.1.4) Generating relationship TG_RZlocusTag TG_oldRZlocusTag
grep  -B1 -wf  <(cut -f1 $output'p1_tmp' | sort | uniq) $annotation'GCF_000268285.2_RPHCH2410v2_genomic.gbff' | sed 's/                     //' | sed 's/--//'| sed '/^$/d' > $output'temporalRelationship.tmp'
python $bin'python/TG-LocusTag-parsey.py' && rm $output'temporalRelationship.tmp'


# 1.1.5) Generating temporary file TF_RegunlobdbID, TG_RegunlobdbID, TG_name, TFname_TGRegunlondbId
# Double check the grep -wf <(cut -f4 "$output/p1_tmp") <(grep -v "#" $networkInfo'network_tf_gene.txt' | cut -f1,2,3,4 | sort | uniq) | cut -f2 | sort | uniq | wc == 16; if not check synonyms
grep -wf <(cut -f9 $OrthologousTFInfo'Orthologous_table.tsv' | sort | uniq) <(paste <(grep -wf <(cut -f4 "$output/p1_tmp") <(grep -v "#" $networkInfo'network_tf_gene.txt' | cut -f1,2,3,4 | sort | uniq) | cut -f1,3,4) <(grep -wf <(cut -f4 "$output/p1_tmp") <(grep -v "#" $networkInfo'network_tf_gene.txt' | cut -f1,2,3,4 | sort | uniq) | cut -f2,3 | sed 's/\t/-/g')) > $output'p5_tmp'


# 2) Assembling TF-TGfirstGene.tsv
Rscript --vanilla $bin'R/oTF-TG-NetworkGenerator.R' $output'p1_tmp' $output'p2_tmp' $output'p3_tmp' $output'p4_tmp' $output'p5_tmp' $output'oTF-TGfirstGene-TRN-bgM1.tsv' && cd $output && rm *tmp

# 3) Expanding the TRN to regulons
python3 $bin'python/extending_trn_to_operons.py'

#########################################################################################################
#################################### Plotting the TRN (Main figure 4) ###################################
#########################################################################################################
Rscript --vanilla $bin'R/section3-MainFigure4.R'  $tables'MainTab1.tmp' $png

#########################################################################################################
#################################### Computing table 1 (Main table 1) ###################################
#########################################################################################################

Rscript --vanilla $bin'R/section3-MainTable1.R' $output'oTF-TGs-TRN-bgM1.tsv' $OrthologousTFInfo'TF_Orthologous_table.tsv' $OrthologousTFInfo'ECaaq_oTFs_motifs_info.tsv' $OrthologousTFInfo'RZaaq_oTFs_motifs_info.tsv' $tables'MainTab1.tsv' && rm $tables'MainTab1.tmp'