# Name:
#  BBH_TFs.sh
# Author:
#  Hernandez-Benitez Ericka Montserrat
# Version
#  v1.7
# 

# Description
#	The script performs the following analyses:
#		1) BBH between E. coli and R phaseoli genomes with blastall program
#		2) Identification of orthologus genes basen on BBH with identity >= 30, bit score >= 50 and e-value <= 0.001
#		3) Gene annotation of the orthologus genes (i.e. locus tag, common name, protein ID, among others)
#		5) Identification of orthologus TFs
# 		6) Plots of BLAT quality (i.e coverage, identity, e-value and qlen)
# sh BBH_TFs.sh


# Additional notes:
# 	The script contemplates the following structure in directories and thus
#	the script will NOT generate them.

#.
#├── annotation
#│   ├── merge_annotation.tsv
#│   └── TFs_coli_v2.txt
#├── bin
#│   ├── bash
#│   │   └── BBH_TFs.sh
#│   ├── python
#│   │   ├── motifs_seq_relation_v2-Parser.py
#│   │   └── RZ_TFs_annotation_parser.py
#│   └── R
#│       ├── merge_coverage_BLAST.R
#│       ├── merge_Etables_p1_v3.R
#│       ├── merge_Etables_p2_v4.R
#│       ├── merge_Etables_p3_v5.R
#│       ├── merge_motif_seq_v2.R
#│       ├── results_section1.R
#│       ├── section1-Figures-main.R
#│       └── section1-Figures-supp.R
#├── genomesInfo
#│   ├── EC
#│   │   ├── GCF_000005845.2_ASM584v2_genomic.fna
#│   │   └── GCF_000005845.2_ASM584v2_protein.faa
#│   └── RZ
#│       ├── GCA_000268285.2_RPHCH2410v2_genomic.gtf
#│       ├── GCF_000268285.2_RPHCH2410v2_cds_from_genomic.fna
#│       ├── GCF_000268285.2_RPHCH2410v2_genomic.fna
#│       ├── GCF_000268285.2_RPHCH2410v2_genomic.fna.fai
#│       ├── GCF_000268285.2_RPHCH2410v2_genomic.gbff
#│       ├── GCF_000268285.2_RPHCH2410v2_genomic.gff
#│       └── GCF_000268285.2_RPHCH2410v2_protein.faa
#├── metabolitesInfo
#│   └── TFs_metabolites.txt
#└── motifsInfo
    #└── motifs_seq_relation.tsv
#

#########################################################################################################
######################################### Working directory #############################################
#########################################################################################################

cd /home/emhernan/1_BBH_TFs/

#########################################################################################################
##################################### Saving all work directories #####################################
#########################################################################################################

indir=/home/emhernan/1_BBH_TFs/
bin=/home/emhernan/1_BBH_TFs/bin/ # must be created with all files
genomesInfo=/home/emhernan/1_BBH_TFs/genomesInfo/ #must be created with genome files
ECdb=/home/emhernan/1_BBH_TFs/ECdb/
RZdb=/home/emhernan/1_BBH_TFs/RZdb/
BLASTresults=/home/emhernan/1_BBH_TFs/BLASTresults/
ortologousInfo=/home/emhernan/1_BBH_TFs/ortologousInfo/
motifsInfo=/home/emhernan/1_BBH_TFs/motifsInfo/ # must be created with the Ecocyc motifs file
annotation=/home/emhernan/1_BBH_TFs/annotation/ # must be created with merge_annotation.tsv and TFs_coli_v2.txt files
tables=/home/emhernan/1_BBH_TFs/tables/ 
metabolitesInfo=/home/emhernan/1_BBH_TFs/metabolitesInfo/ # must be created with TFs_metabolites.txt
png=/home/emhernan/1_BBH_TFs/png/


#########################################################################################################
############################################### BLAST BBH ###############################################
#########################################################################################################

cd $indir
mkdir $ECdb
mkdir $RZdb
mkdir $BLASTresults

# 1) Generate a proper fasta header for ECdb indexing
perl -pe 'if(/^>((NP|YP)_\d+\.\d*)(.+)/){ $c++; s/>/>gnl|ECaadb|$c|$1/ }' $genomesInfo'EC/GCF_000005845.2_ASM584v2_protein.faa' > $ECdb'GCF_000005845.2_ASM584v2_protein.faaed'

# 2) Run formatdb, generating an indexed DB for EC 
formatdb -i $ECdb'GCF_000005845.2_ASM584v2_protein.faaed' -p T -o T

# 3) Generate a proper fasta header for RZdb indexing
perl -pe 'if( /^>(WP_\d+\.\d*)(.+)/ ){ $c++; s/>/>gnl|RZaadb|$c|$1/ }' $genomesInfo'RZ/GCF_000268285.2_RPHCH2410v2_protein.faa' > $RZdb'GCF_000268285.2_RPHCH2410v2_protein.faaed'


# 4) Run formatdb, generating an indexed DB for RZ
formatdb -i $RZdb'GCF_000268285.2_RPHCH2410v2_protein.faaed' -p T -o T


# 5) Run blastp (blastall -p blastp), and get only the best hit (-b 1) in tabular format (-m 8)
# 5.1) RZ as query
blastall -p blastp -i $RZdb'GCF_000268285.2_RPHCH2410v2_protein.faaed' -a 6 -d $ECdb'GCF_000005845.2_ASM584v2_protein.faaed' -b 1 -m 8 -e 0.001 >  $BLASTresults'RZaaq_ECaadb_blastP_b1_m8.tab'
# 5.2) EC as query
blastall -p blastp -i $ECdb'GCF_000005845.2_ASM584v2_protein.faaed' -a 6 -d $RZdb'GCF_000268285.2_RPHCH2410v2_protein.faaed' -b 1 -m 8 -e 0.001 > $BLASTresults'ECaaq_RZaadb_blastP_b1_m8.tab'
# Selenocysteine (U) at position 140 replaced by X

#########################################################################################################
################################ Generating BLAST results with coverage #################################
#########################################################################################################

# 1) Generate length sequences for Escherichia coli
paste <(grep ">" $ECdb'GCF_000005845.2_ASM584v2_protein.faaed'| cut -f1 -d " " | sed 's/>//g') <( perl -nae 'if(/(^[GPAVLIMCFYWHKRQNEDSTUO]*)\n/) {print "$1";} else{ print "\n";}' $ECdb'GCF_000005845.2_ASM584v2_protein.faaed' | awk '{print length}' |  sed -n '2,4287p') > $BLASTresults'ECgene_len.tsv'

# 2) Generate length sequences for Rhizobium phaseoli
paste <(grep ">" $RZdb'GCF_000268285.2_RPHCH2410v2_protein.faaed' | cut -f1 -d " " | sed 's/>//g') <(perl -nae 'if(/(^[GPAVLIMCFYWHKRQNEDSTUO]*)\n/) {print "$1";} else{ print "\n";}' $RZdb'GCF_000268285.2_RPHCH2410v2_protein.faaed' | awk '{print length}'| sed -n '2,6198p') > $BLASTresults'RZgene_len.tsv'

# 3) Merge BLAST v1 with genLen 
cd $BLASTresults
Rscript --vanilla $bin'R/merge_coverage_BLAST.R' $BLASTresults'ECaaq_RZaadb_blastP_b1_m8.tab' $BLASTresults'ECgene_len.tsv' V1 V1 $BLASTresults'ECaaq_RZaadb_blastP_b1_m8_v2.tab'
Rscript --vanilla $bin'R/merge_coverage_BLAST.R' $BLASTresults'RZaaq_ECaadb_blastP_b1_m8.tab' $BLASTresults'RZgene_len.tsv' V1 V1 $BLASTresults'RZaaq_ECaadb_blastP_b1_m8_v2.tab'

#########################################################################################################
############################### Generating gene list with motif sequences ###############################
#########################################################################################################

# 1) Generating motifs_seq_relation_v2.tsv table: geneName	locusTag	motifDesc	start	end	motifSeq
python $bin'python/motifs_seq_relation_v2-Parser.py'

# 2) FIltering the motifs we are interested 
grep -wE "Ca-Binding-Region|Conserved-Region|Catalytic-Domain|DNA-Binding-Region|Intramembrane-Region|Nucleotide-Phosphate-Binding-Region|Protein-Structure-Region|Alpha-Helix-Region|Beta-Strand-Region|Coiled-Coil-Region|Transmembrane-Region|Zn-Finger-Region" $motifsInfo'motifs_seq_relation_v2.tsv' > ed && mv ed $motifsInfo'motifs_seq_relation_v2.tsv'

# 3) Generating gene_aa_seq.tsv table: protein_ID	sequence
paste <(grep ">" $ECdb'GCF_000005845.2_ASM584v2_protein.faaed'| cut -f1 -d " " | sed 's/>//g' | cut -f4 -d "|") <(perl -nae 'if(/(^[GPAVLIMCFYWHKRQNEDSTUO]*)\n/) {print "$1";} else{ print "\n";}' $ECdb'GCF_000005845.2_ASM584v2_protein.faaed'| sed -n '2,4286p') > $motifsInfo'gene_aa_seq.tsv'

# 4) Generating a prot_lt_tmp temporary table: protein_ID	LT --> Preguntar a Dani si se debe poner el codigo de la generacion de la anotacion
cut -f1,2 $annotation'merge_annotation.tsv' > $motifsInfo'prot_lt_tmp'

# 5) Merging all the previous table to generate a tsv file named motifs_seq_relation_v3.tsv: locusTag	motifDesc	mSS	mSE	motifsSeq
Rscript --vanilla $bin'R/merge_motif_seq_v2.R' $motifsInfo $motifsInfo'motifs_seq_relation_v3.tsv' gene_aa_seq.tsv motifs_seq_relation_v2.tsv prot_lt_tmp && rm $motifsInfo'prot_lt_tmp'

#########################################################################################################
######################################### Filtering Orthologous #########################################
#########################################################################################################

mkdir $ortologousInfo

# 1) Retriving Orthologous BLAST ID sequences: with a 0.001 e-value cut-off, bit score >= 50, and identity >= 30

ECqRZdb_ort=$(paste <(awk '{if ($11 <= 0.00000023 && $3 >=30 && $12 >=50 && $14 >=80) {print }}' $BLASTresults'ECaaq_RZaadb_blastP_b1_m8_v2.tab' | cut -f1 | sed -e 's/|[YNP]*_[0-9]*.[0-9]*//g' | cut -f1) <(awk '{if ($11 <= 0.00000023 && $3 >=30 && $12 >=50 && $14 >=80) {print }}' $BLASTresults'ECaaq_RZaadb_blastP_b1_m8_v2.tab' | cut -f2))

RZqECdb_ort=$(paste <(awk '{if ($11 <= 0.00000016 && $3 >=30 && $12 >=50 && $14 >=80) {print }}' $BLASTresults'RZaaq_ECaadb_blastP_b1_m8_v2.tab' | cut -f1 | sed 's/|WP_[0-9]*.[0-9]*//g') <(awk '{if ($11 <= 0.00000016 && $3 >=30 && $12 >=50 && $14 >=80) {print }}' $BLASTresults'RZaaq_ECaadb_blastP_b1_m8_v2.tab' | cut -f2))

intersec_ort=$(grep -wf <( awk '{if ($11 <= 0.00000023 && $3 >=30 && $12 >=50 && $14 >=80) {print }}' $BLASTresults'ECaaq_RZaadb_blastP_b1_m8_v2.tab'| cut -f2) <( awk '{if ($11 <= 0.00000016 && $3 >=30 && $12 >=50 && $14 >=80) {print }}' $BLASTresults'RZaaq_ECaadb_blastP_b1_m8_v2.tab' | cut -f1 | sed 's/|WP_[0-9]*.[0-9]*//g'))

grep -wf <(echo "$intersec_ort") <(echo "$ECqRZdb_ort") | wc #1217

# 1266 #975 Orthologous genes
grep -wf <(grep -wf <(echo "$intersec_ort") <(echo "$RZqECdb_ort") | perl -nae 'print "$F[1]\t$F[0]\n"') <(echo "$ECqRZdb_ort") | sort |uniq  > $ortologousInfo'allOrthologousID.tsv'

# Alternative
intersec2_ort=$(grep -wf <(cut -f1 $BLASTresults'RZaaq_ECaadb_blastP_b1_m8_v2.tab' | sed 's/|WP_[0-9]*.[0-9]*//g') <(cut -f2  $BLASTresults'ECaaq_RZaadb_blastP_b1_m8_v2.tab')) # 4148 #2443
grep -wf <(grep -wf <(echo "$intersec2_ort") <(echo "$ECqRZdb_ort") | perl -nae 'print "$F[1]\t$F[0]\n"') <(echo "$RZqECdb_ort") | sort |uniq | wc # 1266 # 975
#########################################################################################################
####################################### Retrieving BLAST sequences ######################################
#########################################################################################################

mkdir $tables

# 1) Orthologous table

# 1.1) Generating Orthologous_ECaaq_RZaadb_blastP_b1_m8.tab
grep -wf <(cut -f2 $ortologousInfo'allOrthologousID.tsv') <(grep -wf <(cut -f1 $ortologousInfo'allOrthologousID.tsv') $BLASTresults'ECaaq_RZaadb_blastP_b1_m8_v2.tab' | awk '{if ($11 <= 0.00000016 && $3 >=30 && $12 >=50 && $14 >=80) {print }}') > $BLASTresults'Orthologous_ECaaq_RZaadb_blastP_b1_m8.tab'
grep -wf <(cut -f2  $BLASTresults'Orthologous_ECaaq_RZaadb_blastP_b1_m8.tab') <(grep -wf <(cut -f1 $BLASTresults'Orthologous_ECaaq_RZaadb_blastP_b1_m8.tab' | cut -f1,2,3 -d "|") $BLASTresults'RZaaq_ECaadb_blastP_b1_m8_v2.tab' | awk '{if ($11 <= 0.00000023 && $3 >=30 && $12 >=50 && $14 >=80) {print }}') > $BLASTresults'Orthologous_RZaaq_ECaadb_blastP_b1_m8.tab'
cat <(head -n1 $BLASTresults'ECaaq_RZaadb_blastP_b1_m8_v2.tab') $BLASTresults'Orthologous_ECaaq_RZaadb_blastP_b1_m8.tab' > ed && mv ed $BLASTresults'Orthologous_ECaaq_RZaadb_blastP_b1_m8.tab'
cat <(head -n1 $BLASTresults'RZaaq_ECaadb_blastP_b1_m8_v2.tab') $BLASTresults'Orthologous_RZaaq_ECaadb_blastP_b1_m8.tab'  > ed && mv ed $BLASTresults'Orthologous_RZaaq_ECaadb_blastP_b1_m8.tab'


# 1.2)  Generating Orthologous_ECaaq_RZaad_BLASTid_protein.tsv table: EC_BLASID	RZ_BLASTID	EC_proteinID
paste <(cut -f1,2 $BLASTresults'Orthologous_ECaaq_RZaadb_blastP_b1_m8.tab') <(cut -f1 $BLASTresults'Orthologous_ECaaq_RZaadb_blastP_b1_m8.tab'| cut -f4 -d "|") | grep -v "qName" | sort | uniq > $tables'Orthologous_ECaaq_RZaad_BLASTid_protein.tsv'

# 1.2)  Generating Orthologous_RZaad_ECaaq_BLASTid_protein.tsv table: RZ_BLASTID	EC_BLASID_inc	RZ_proteinID
paste <(paste <(cut -f1 $BLASTresults'Orthologous_RZaaq_ECaadb_blastP_b1_m8.tab') <(cut -f1 $BLASTresults'Orthologous_RZaaq_ECaadb_blastP_b1_m8.tab'| cut -f4 -d "|")) <(cut -f1 $BLASTresults'Orthologous_RZaaq_ECaadb_blastP_b1_m8.tab' | cut -f1,2,3 -d "|") | grep -v "qName" | sort | uniq  > $tables'Orthologous_RZaad_ECaaq_BLASTid_protein.tsv'

# 1.3)  Generating RZ_functional_annotation.tsv table: RZ_locusTag	RZ_product	RZ_proteinID (annotation for all RZ genes)
cat $genomesInfo'RZ/GCF_000268285.2_RPHCH2410v2_genomic.gbff' | grep "/" | sed -r 's/(.*)\/(.*)=(.*)/\2=\3/' |  grep -wE "locus_tag|product|protein_id" |perl -nae 'if(/(.*)\="(RPH.*)"/){print"\n$2\t"} else{if(/(.*)\="(.*)\"?/){print"$2\t"}}' | sed 's/"//g' | perl -nae 'if(/(.*)\t(.*)\t(.*)/){print "$1\t$2\t$3\n"}' | cut -f1,2,3 > $annotation'RZ_functional_annotation.tsv'

# 1.4)  Generating Orthologous_RZ_annotation.tsv table: RZ_locusTag	RZ_product	RZ_proteinID (annotation only for orthologous)
grep -wf <(cut -f1 $BLASTresults'Orthologous_RZaaq_ECaadb_blastP_b1_m8.tab'  | cut -f4 -d "|"| sort | uniq) $annotation'RZ_functional_annotation.tsv' | sed s/\'//g > $tables'Orthologous_RZ_annotation.tsv'

# 1.5)  Generating BLAST_ECq_tmp temporary table: EC_BLASID	qSS	qSE
cut -f1,7,8 $BLASTresults'Orthologous_ECaaq_RZaadb_blastP_b1_m8.tab' > $tables'BLAST_ECq_tmp'

# 1.6)	Merging table part_1: ECBLAST_ID	RZBLAST_ID	EC_locusTag	NCBI_name	Regulondb_name	Abasy_name	Ecocyc_name	Synonyms	RegulondbID	EC_proteinID	RZ_locusTag	RZ_proteinID	qSS	qSE	BLASTseq	EC_product	RZ_Product
Rscript --vanilla $bin'/R/merge_Etables_p1_v3.R' $tables $tables'Orthologous_table_p1.tsv' Orthologous_ECaaq_RZaad_BLASTid_protein.tsv /home/emhernan/1_BBH_TFs/annotation/merge_annotation.tsv Orthologous_RZaad_ECaaq_BLASTid_protein.tsv Orthologous_RZ_annotation.tsv $motifsInfo'gene_aa_seq.tsv' BLAST_ECq_tmp && rm $tables'BLAST_ECq_tmp'

# 1.7)	Merging table part_2: ECBLAST_ID	RZBLAST_ID	EC_locusTag	NCBI_name	Regulondb_name	Abasy_name	Ecocyc_name	Synonyms	RegulondbID	EC_proteinID	RZ_locusTag	RZ_proteinID	qSS	qSE	BLASTseq	EC_product	RZ_Product	motifDesc	mSS	mSE	motifsSeq
# Note: Motif_conservation will be done in 2_MotifAnalysis

Rscript --vanilla $bin'/R/merge_Etables_p2_v4.R' $motifsInfo'motifs_seq_relation_v3.tsv' $tables'Orthologous_table_p1.tsv' $tables'Orthologous_table_p2.tsv' && rm $tables'Orthologous_table_p1.tsv'

# 1.8)	Generating Metabolites_relation_tmp temporary table: EC_locusTag	Metabolites
grep -v "#" $metabolitesInfo'TFs_metabolites.txt' | cut -f2,4 | sed s/\'//g > $tables'Metabolites_relation_tmp'

# 1.9) Generating old RZ locus tag, new RZ locus tag for all RZ genes
cat $genomesInfo'RZ/GCF_000268285.2_RPHCH2410v2_genomic.gbff' | grep "/" | sed -r 's/(.*)\/(.*)=(.*)/\2=\3/' |  grep -wE "locus_tag|old_locus_tag|product" |perl -nae 'if(/(product)\="(.*)\"?/){print"$2\n"} else{if(/(.*)\="(.*)\"?/){print"$2\t"}}' | sed 's/"//g' | cut -f1,2 > $annotation'/RZ_functional_annotation_v2.tsv'

# 1.10)	Merging table part_3: ECBLAST_ID	RZBLAST_ID	EC_locusTag	NCBI_name	Regulondb_name	Abasy_name	Ecocyc_name	Synonyms	RegulondbID	EC_proteinID	RZ_locusTag	RZ_proteinID	qSS	qSE	BLASTseq	EC_product	RZ_Product	motifDesc	mSS	mSE	motifsSeq	effector_name 
Rscript --vanilla $bin'R/merge_Etables_p3_v5.R' $tables'Metabolites_relation_tmp' $annotation'RZ_functional_annotation_v2.tsv' $tables'Orthologous_table_p2.tsv' $tables'Orthologous_table.tsv' && cd $tables && rm  Orthologous_table_p2.tsv && rm *tmp

# 2) TFs table
cat <(head -n1 $tables'Orthologous_table.tsv') <(grep -wf $annotation'TFs_coli_v2.txt' <(cat $tables'Orthologous_table.tsv')) > $tables'TF_Orthologous_table.tsv'

# 3) TFs table from BLAST

cat <(head -n1 $BLASTresults'Orthologous_ECaaq_RZaadb_blastP_b1_m8.tab') <(grep -wf <(cut -f1 $tables'TF_Orthologous_table.tsv' | sort | uniq) $BLASTresults'Orthologous_ECaaq_RZaadb_blastP_b1_m8.tab') > $BLASTresults'TFs_ECaaq_RZaadb_blastP_b1_m8.tab'
cat <(head -n1 $BLASTresults'Orthologous_RZaaq_ECaadb_blastP_b1_m8.tab') <(grep -wf <(cut -f2 $tables'TF_Orthologous_table.tsv' | sort | uniq) $BLASTresults'Orthologous_RZaaq_ECaadb_blastP_b1_m8.tab') > $BLASTresults'TFs_RZaaq_ECaadb_blastP_b1_m8.tab'


#########################################################################################################
######################################### Plotting Main Figure 1  #######################################
#########################################################################################################
mkdir $png

Rscript --vanilla $bin'R/section1-Figures-main.R' $BLASTresults $png
Rscript --vanilla $bin'R/section1-Figures-supp.R' $BLASTresults $png

#########################################################################################################
######################## Checking the number of TFs previously annotated as TFs #########################
#########################################################################################################

cd $tables

# 1) Search the annotation using the protein ID. The annotitation come from GCF_000268285.2_RPHCH2410v2_genomic.gbff file

grep -B3 -wf <(grep -wf <(cut -f2 $tables'TF_Orthologous_table.tsv' | sort | uniq | grep -v "RZBLAST_ID" | cut -f4 -d "|") $genomesInfo'RZ/GCF_000268285.2_RPHCH2410v2_genomic.gbff' | grep -w "protein_id") $genomesInfo'RZ/GCF_000268285.2_RPHCH2410v2_genomic.gbff' | grep -v "transl_table" | grep -v "codon_start"> $tables'TF_Orthologous_RZannotation_table.tsv'
python $bin'python/RZ_TFs_annotation_parser.py'

grep -wE "regulator|activator|repressor|regulation protein| regulatory protein" $tables"TF_Orthologous_RZannotation_table_v2.tsv" | wc 
# Remember 57 oTFs. 