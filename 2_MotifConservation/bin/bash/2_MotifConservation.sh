#########################################################################################################
######################################### Working directory #############################################
#########################################################################################################

cd /home/emhernan/2_MotifConservation


#########################################################################################################
##################################### Saving all work directories #####################################
#########################################################################################################

indir=/home/emhernan/2_MotifConservation
bin=/home/emhernan/2_MotifConservation/bin/ # must be created with all files
genomesInfo=/home/emhernan/2_MotifConservation/genomesInfo/ #must be created with genome files
BLASTresultsSeq=/home/emhernan/2_MotifConservation/BLASTresultsSeq/
ECdb=/home/emhernan/2_MotifConservation/ECdb/
RZdb=/home/emhernan/2_MotifConservation/RZdb/
formated=/home/emhernan/2_MotifConservation/formated/
motifsInfo=/home/emhernan/2_MotifConservation/motifsInfo/
png=/home/emhernan/2_MotifConservation/png/



#########################################################################################################
############################################### BLAST BBH ###############################################
#########################################################################################################
cd $indir
mkdir $ECdb
mkdir $RZdb
mkdir $BLASTresultsSeq

# 1) Generate a proper fasta header for ECdb indexing
perl -pe 'if(/^>((NP|YP)_\d+\.\d*)(.+)/){ $c++; s/>/>gnl|ECaadb|$c|$1/ }' $genomesInfo'EC/GCF_000005845.2_ASM584v2_protein.faa' > $ECdb'GCF_000005845.2_ASM584v2_protein.faaed'

# 2) Run formatdb, generating an indexed DB for EC 
formatdb -i $ECdb'GCF_000005845.2_ASM584v2_protein.faaed' -p T -o T

# 3) Generate a proper fasta header for RZdb indexing
perl -pe 'if( /^>(WP_\d+\.\d*)(.+)/ ){ $c++; s/>/>gnl|RZaadb|$c|$1/ }' $genomesInfo'RZ/GCF_000268285.2_RPHCH2410v2_protein.faa' > $RZdb'GCF_000268285.2_RPHCH2410v2_protein.faaed'


# 4) Run formatdb, generating an indexed DB for RZ
formatdb -i $RZdb'GCF_000268285.2_RPHCH2410v2_protein.faaed' -p T -o T


# 5) Run blastp (blastall -p blastp), and get only the best hit (-b 1) with pairwise format (m0)

blastall -p blastp -i $RZdb'GCF_000268285.2_RPHCH2410v2_protein.faaed' -a 6 -d $ECdb'GCF_000005845.2_ASM584v2_protein.faaed' -b 1 -m 0 >  $BLASTresultsSeq'RZaaq_ECaadb_blastP_b1_m0_seq.txt'
# 5.2) EC as query
blastall -p blastp -i $ECdb'GCF_000005845.2_ASM584v2_protein.faaed' -a 6 -d $RZdb'GCF_000268285.2_RPHCH2410v2_protein.faaed' -b 1 -m 0 > $BLASTresultsSeq'ECaaq_RZaadb_blastP_b1_m0_seq.txt'
# Selenocysteine (U) at position 140 replaced by X


#########################################################################################################
###################################### Extracting BLAST sequences  ######################################
#########################################################################################################

mkdir $formated

# 1) Extracting sequences from ECaaq_RZaadb_blastP_b1_m0_seq.txt to create ECaaq_RZaadb_blastP_seq_formated.tsv and
#  ECaaq_RZaadb_blastP_seq_formated_v2.tsv files with the following columns: start	end	seq	seqIdent (with and without gaps)

# 1.1) Generate ECaaq_RZaadb_blastP_seq.tab file: start, sequence, end (the sequences per ID are not in one line)
grep -E -A1 "Query:|Score = " $BLASTresultsSeq'ECaaq_RZaadb_blastP_b1_m0_seq.txt' | sed 's/Query: //' | sed 's/Query: //' | grep -v "Identities" | sed 's/ Score = /Score\t/' | perl -nae 'if (/^(Score)\t.*/) {print "$1\t\t\n"}; if(/^(\d+)\s+(.+)\s+(\d+)/){print "$1\t$2\t$3\n"}; ; if(/^\s+(.*)/) {print "\t$1\t\n"}' | tail -n +2  > $formated'ECaaq_RZaadb_blastP_seq.tab'

# 1.2) Generate ECaaqID_temp.tsv: Sequence identifiers repited n times. N = number of best aligments per idenfier
grep -E "Query= | Score ="  $BLASTresultsSeq'ECaaq_RZaadb_blastP_b1_m0_seq.txt' | sed 's/Query= //' | sed 's/ Score = /Score\t\t/g' | sed 's/ /\t/g' | cut -f1 > $formated'ECaaqID_temp.tsv'

# 1.3) To create ECaaq-subjID.tsv file: IDs from the subject (IDs from R. phaseoli sequences)
conda activate
python $bin'python/ECaaq-BLASTseqParser-subjID.py'

# 1.4) To create ECaaq-queID.tsv file: IDs from the query (IDs from Ecoli sequences)
python $bin'python/ECaaq-BLASTseqParser-queID.py'


# 1.5) To create ECaaq_RZaadb_blastP_seq_formated.tsv: start	end	seq	seqIdent (with gaps)
python $bin'python/ECaaq-BLASTseqParser.py'

# 1.6) To create ECaaq_RZaadb_blastP_seq_formated_v2.tsv: start	end	seq	seqIdent (withOutGaps)
python $bin'python/ECaaq-MiddleDash-parser.py'
cd $formated


# 2) Extracting sequences from RZaaq_ECaadb_blastP_b1_m3_seq.txt to create RZaaq_ECaadb_blastP_seq_formated.tsv and
#  RZaaq_ECaadb_blastP_seq_formated_v2.tsv files with the following columns: start	end	seq	seqIdent (with and without gaps)


# 2.1) Generate ECaaqID_temp.tsv file:  start, sequence, end (the sequences per ID are not in one line)
grep -E -B1 "Sbjct:|Score = " $BLASTresultsSeq'RZaaq_ECaadb_blastP_b1_m0_seq.txt' | sed 's/Sbjct: //' | grep -v "Identities" | sed 's/ Score = /Score\t/' | perl -nae 'if (/^(Score)\t.*/) {print "$1\t\t\n"}; if(/^(\d+)\s+(.+)\s+(\d+)/){print "$1\t$2\t$3\n"}; ; if(/^\s+(.*)/) {print "\t$1\t\n"}' | tail -n +3  > $formated'RZaaq_ECaadb_blastP_seq.tab'

# 2.2) Generate RZaaqID_temp.tsv: Sequence identifiers repited n times. N = number of best aligments per idenfier
grep -E "Query= | Score ="  $BLASTresultsSeq'RZaaq_ECaadb_blastP_b1_m0_seq.txt' | sed 's/Query= //' | sed 's/ Score = /Score\t\t/g' | sed 's/ /\t/g' | cut -f1 > $formated'RZaaqID_temp.tsv'


# 2.3) To create RZaaq-queID.tsv: IDs from the subject (IDs from Ecoli sequences)
python $bin'python/RZaaq-BLASTseqParser-queID.py'

# 2.4) To create RZaaq_ECaadb_blastP_seq_formated.tsv: start	end	seq	seqIdent (with gaps)
python $bin'python/RZaaq-BLASTseqParser.py'

# 2.5) To create RZaaq_ECaadb_blastP_seq_formated_v2.tsv: start	end	seq	seqIdent (withOutGaps)
python $bin'python/RZaaq-MiddleDash-parser.py'
cd $formated


#########################################################################################################
################################### Computing conservation of motifs  ###################################
#########################################################################################################

mkdir $motifsInfo

# 1) Compute identity and coverage of both BLAST for ALL motifs 
# Obtain ECaaq_oTFs_motifs_info.tsv and RZaaq_oTFs_motifs_info.tsv files
sed -i '2d' ECaaq_RZaadb_blastP_seq_formated.tsv
cp /home/emhernan/1_BBH_TFs/tables/TF_Orthologous_table.tsv $formated 

Rscript --vanilla $bin'R/AllMotifs-processing.R' $formated $motifsInfo 'ECaaq-queID.tsv' 'ECaaq_RZaadb_blastP_seq_formated.tsv' 'ECaaq_RZaadb_blastP_seq_formated_v2.tsv' 'ECBLAST_ID' 'ECaaq_oTFs_motifs_info.tsv'
Rscript --vanilla $bin'R/AllMotifs-processing.R' $formated $motifsInfo 'RZaaq-queID.tsv' 'RZaaq_ECaadb_blastP_seq_formated.tsv' 'RZaaq_ECaadb_blastP_seq_formated_v2.tsv' 'RZBLAST_ID' 'RZaaq_oTFs_motifs_info.tsv'

# 2) Perform conservation analysis: motif conservation if identity >= 30
Rscript --vanilla $bin'R/motifConservationAnalysis.R' $motifsInfo

# 3) Ploting pie chart (Supplementary Figure 3)
mkdir $png
Rscript --vanilla $bin'R/section2-piechartFigure.R' $motifsInfo'motifsConservationIdent30.tsv' $png

# 4) Filtering DNA-binding motifs results from the conservation analysis

grep "DNA-Binding-Region" $motifsInfo'ECaaq_oTFs_motifs_info.tsv' > $motifsInfo'ECaaq_oTFs_DNAbinding_motifs_info.tsv'
grep "DNA-Binding-Region" $motifsInfo'RZaaq_oTFs_motifs_info.tsv' > $motifsInfo'RZaaq_oTFs_DNAbinding_motifs_info.tsv'

# 5) Plotting identity vs coverage of DNA binding motifs from oTFs (Main Figure 2)
Rscript --vanilla $bin'R/section2-MainFigure2.R' $motifsInfo $png ECaaq_oTFs_DNAbinding_motifs_info.tsv RZaaq_oTFs_DNAbinding_motifs_info.tsv mako MainF2.png TRUE

