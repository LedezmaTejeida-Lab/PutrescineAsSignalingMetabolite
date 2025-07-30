docs=/space24/PGC/emhernan/4_1_puuRPSSM_evaluation/docs/
bin=/space24/PGC/emhernan/4_1_puuRPSSM_evaluation/bin/
matrixes_regulonDB=/space24/PGC/emhernan/4_1_puuRPSSM_evaluation/matrixes_regulonDB/
rsat=/space24/PGC/emhernan/4_1_puuRPSSM_evaluation/rsat/
sequences=/space24/PGC/emhernan/4_1_puuRPSSM_evaluation/rsat/sequences/
matrix_scan=/space24/PGC/emhernan/4_1_puuRPSSM_evaluation/rsat/matrix_scan/
matrix_scan_regulonDB=/space24/PGC/emhernan/4_1_puuRPSSM_evaluation/rsat/matrix_scan/matrix_scan_regulonDB/
bgfilePath=/space23/rsat/rsat/public_html/data/genomes/Escherichia_coli_str._K-12_substr._MG1655_GCF_000005845.2_ASM584v2/oligo-frequencies/2nt_upstream-noorf_Escherichia_coli_str._K-12_substr._MG1655_GCF_000005845.2_ASM584v2-ovlp-1str.freq.gz
output_PSSM_RegulonDB=/space24/PGC/emhernan/4_1_puuRPSSM_evaluation/output_PSSM_RegulonDB/

mkdir -p /space24/PGC/emhernan/4_1_puuRPSSM_evaluation/
mkdir -p $docs
mkdir -p $bin
mkdir -p $bin'R'
mkdir -p $bin'python'
mkdir -p $rsat
mkdir -p $sequences
mkdir -p $matrixes_regulonDB
mkdir -p $matrix_scan
mkdir -p $output_PSSM_RegulonDB
mkdir -p $matrix_scan_regulonDB

# Download GeneProductAllIdentifiersSet.txt, TFSet.txt and PSSM_v4_0 files from https://testregulondb.ccg.unam.mx/menu/download/datasets/index.jsp
# RegulonDB 11.2, 05/22/2023.
# These files were used to map the locus-tag of the TF name of the PSSM_v4_0 files

# Download TF-RISet.tsv GeneProductSet.tsv and OperonSet.tsv  from RegulonDB Version: 13.5.0 01-28-2025
# Release Date: 2025-01-20
# These files were used to map the locus-tag of the TG name of the TF-RISet file


# 0) Compute the number of  oTFs with coonserved DBM #49
cp -r /space24/PGC/emhernan/docs/4_1_puuRPSSM_evaluation/* $docs
cp /space24/PGC/emhernan/bin/R/4_1_puuRPSSM_evaluation/* $bin'R'
cp /space24/PGC/emhernan/bin/python/4_1_puuRPSSM_evaluation/* $bin'python'

###################################################################################################################
####################################### PSSM v4 RegulonDB analysis ################################################
###################################################################################################################

# 1) Computing the number of PSSMs for v4
grep "Transcription Factor Name:" $docs'PSSM_v4_0/results/PSSM-Dataset-v4.0.txt' | sort | uniq | wc
# 92 PSSM in v4, this version contains weak and strong matrixes mixed. 

# 2) Extrating the information of PuuR PSSM
grep "Transcription Factor Name:" $docs'PSSM_v4_0/results/PSSM-Dataset-v4.0.txt'| cut -f2 -d ":" | sed 's/^[ \t]*//;s/[ \t]*$//' | grep PuuR > $docs'TFname_PSSMv4.tsv'
grep -v "#" $docs'TFSet.txt' | cut -f 1,2,4 | sort  |uniq | sed '/^$/d' | grep PuuR > $docs'TF_info_v4.txt'
grep -v "#" $docs'GeneProductAllIdentifiersSet.txt' | cut -f2,6 | grep -wE "acrR|agaR|araC|arcA|argP|argR|asnC|baeR|basR|caiF|cpxR|cra|crp|csgD|cysB|cytR|ttdR|dcuR|deoR|dnaA|evgA|exuR|fadR|fhlA|fis|flhC|flhD|fnr|fur|gadE|gadW|gadX|galR|galS|gcvA|glpR|glrR|gntR|hipA|hipB|iclR|ihfA|ihfB|iscR|leuO|lexA|lrp|malT|marA|melR|metJ|metR|mlc|mlrA|mntR|modE|mqsA|mraZ|nac|nagC|nanR|narL|narP|nhaR|nrdR|nsrR|glnG|ompR|oxyR|pdhR|phoB|phoP|purR|putA|puuR|qseB|rcdA|rcsA|rcsB|rcsB|relB|rhaS|rob|rstA|rutR|slyA|soxR|soxS|torR|trpR|tyrR|ulaR|uxuR|xylR|ydeO|yjjQ" | perl -nae 'if(/^(\S+)\t.*(b\d+)\D+/){print "$1\t$2\n"}' | grep puuR > $docs'geneannotation.tsv'
cat $docs'PSSM_v4_0/results/PSSM-Dataset_v4.0.cons' | perl -nae 'if(/.*\s+TFName\:\s+(\w+)\,\s+.*TFBSs\:\s+(\d+)\,.*size\:\s+(\d+)/){print "$1\t$2\t$3\n"}' | grep PuuR > $docs'PSSM_info.tmp'
Rscript --vanilla $bin'R/PSSM_info.R' -t $docs'TFname_PSSMv4.tsv' -i $docs'TF_info_v4.txt' -g $docs'geneannotation.tsv' -p $docs'PSSM_info.tmp' -o $docs'PSSM_info.tsv' && rm $docs'TFname_PSSMv4.tsv' && rm $docs'TF_info_v4.txt' && rm $docs'geneannotation.tsv' && rm $docs'PSSM_info.tmp'

# 3) Retrieving upstream sequences from -400 tp +50 -noorf
# To work with a genome for which no mRNA is available, choose ‘CDS’ as feature type and upstream sequences will be retrieved relative to the start codon
# Change organism for Rhizobium_phaseoli_GCF_000268285.2_RPHCH2410v2
rsat retrieve-seq -org Escherichia_coli_str._K-12_substr._MG1655_GCF_000005845.2_ASM584v2  -feattype CDS -type upstream -format fasta -label id,name -from -400 -to +50 -noorf -all > $sequences'retrieve-seq-Escherichia_coli_str._K-12_substr._MG1655_GCF_000005845.2_ASM584v2-10-50.fna'
# 3.1) Convert file to fasta format with Convert-seq
rsat convert-seq  -i $sequences'retrieve-seq-Escherichia_coli_str._K-12_substr._MG1655_GCF_000005845.2_ASM584v2-10-50.fna' -from  fasta -to fasta  -mask non-dna -o $sequences'retrieve-seq-Escherichia_coli_str._K-12_substr._MG1655_GCF_000005845.2_ASM584v2-10-50.fnaed'
# 3.2) Computing the number of bp for all the retrieved sequences 
grep ">" $sequences'retrieve-seq-Escherichia_coli_str._K-12_substr._MG1655_GCF_000005845.2_ASM584v2-10-50.fnaed' | cut -f3 -d ";" | sed 's/size: //'| awk '{s+=$1} END {print s}' #798584 pb

# 4) Filtering upstream sequences to operon leader genes.
# 4.1) Formatig operon file from regulonDB
grep -v "#" $docs'GeneProductSet.tsv' | cut -f2,3 | sort  |uniq | sed 's/[0-9]*)//g' > $docs'gene-bnumber.tsv'
grep -v "#" OperonSet.tsv  | cut -f7 |  grep -v "7)operonGenes" > "$docs/operon_genes.tsv"
Rscript --vanilla $bin'R/operon_formatted_generator.R' -b $docs'gene-bnumber.tsv' -p "$docs/operon_genes.tsv" -o "$docs/eco.ope"

# 4.2) Extracting head of operons and filtering retrieve seq to operon leader genes
awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq=seq $0} END {if (seq) print seq}' $sequences'retrieve-seq-Escherichia_coli_str._K-12_substr._MG1655_GCF_000005845.2_ASM584v2-10-50.fnaed' > $sequences'retrieve-seq-Escherichia_coli_str._K-12_substr._MG1655_GCF_000005845.2_ASM584v2-10-50.fnaed.tmp'
grep -A1 -wf  <(grep -wf <(cut -f1 -d " " "$docs/eco.ope"  | sort | uniq | sed 's/eco-//' | grep -v "None") <(cut -f1 "$sequences/retrieve-seq-Escherichia_coli_str._K-12_substr._MG1655_GCF_000005845.2_ASM584v2-10-50.fnaed.tmp" | sed 's/>//g'))  "$sequences/retrieve-seq-Escherichia_coli_str._K-12_substr._MG1655_GCF_000005845.2_ASM584v2-10-50.fnaed.tmp" | grep -v "\-\-" > $sequences'retrieve-seq-Escherichia_coli_str._K-12_substr._MG1655_GCF_000005845.2_ASM584v2_head_operons.fnaed' && rm $sequences'retrieve-seq-Escherichia_coli_str._K-12_substr._MG1655_GCF_000005845.2_ASM584v2-10-50.fnaed.tmp' &
#pause; background process; wc == 4816

# 4.3) Computing the number of pb that scan matrix will scan
grep ">" $sequences'retrieve-seq-Escherichia_coli_str._K-12_substr._MG1655_GCF_000005845.2_ASM584v2_head_operons.fnaed'  | cut -f3 -d ";" | sed 's/ size: //' | awk '{s+=$1} END {print s}'
grep -v ">" $sequences'retrieve-seq-Escherichia_coli_str._K-12_substr._MG1655_GCF_000005845.2_ASM584v2_head_operons.fnaed' | awk '{print length($0)}' | paste -sd+ | bc
# 614,222

# 4) Formating PuuR PSSM to tab format
grep -v "#" $docs'PSSM_v4_0/results/PSSM-Dataset-v4.0.txt' | perl -nae 'if(/Transcription Factor Name: (\w+)/) {print ";$1\tTranscription Factor\n"} else {if(/([ACG])(\t[\d+\t+]*)/) {print "$1$2\n"} else {if(/(T)(\t[\d+\t+]*)/) {print "$1$2\n//\n"}}}' | grep -A5 PuuR > "$matrixes_regulonDB/PuuR-PSSM-v4.0-ed.txt"

# 5) Running matrix-scan with the threshold established in PSSM-ThreholdsbyTF_v4.0.txt for PuuR

cd $matrixes_regulonDB
for matrix in $(ls "$matrixes_regulonDB")
do
    TF=$(echo "$matrix" | cut -f1 -d "-")
    Threshold=$(grep -w $TF $docs'PSSM_v4_0/results/PSSM-ThreholdsbyTF_v4.0.txt' | cut -f2  | tr -d '[:space:]')
    t="$TF-$Threshold-bgM1.txt"
    rsat matrix-scan -v 1 -matrix_format tab -m "$matrix" -consensus_name -pseudo 1 -decimals 1 -2str -origin end -bgfile $bgfilePath -bg_pseudo 0.01 -return limits -return sites -return pval -uth pval "${Threshold}" -return rank -i  $sequences'retrieve-seq-Escherichia_coli_str._K-12_substr._MG1655_GCF_000005845.2_ASM584v2_head_operons.fnaed' -seq_format fasta -n score > $matrix_scan_regulonDB$t &

done

# 8) Filtering TF-TG interactions for all type of interactions (tf-promoter tf-tu, tf-gene) and the first gene to be transcribed
grep -v "#"  $docs'TF-RISet.tsv' | cut -f2,4,17,20 | sort | uniq | grep -v "riType"  > $docs'TF-TG-alltype-firstGene.tsv'

# 8.1) Mapping first gene of TU name to bnumber
grep -v "#" $docs'GeneProductSet.tsv' | cut -f2,3 | sort  |uniq | sed 's/[0-9]*)//g' > $docs'gene-bnumber.tsv'

# 9) Generating TF-TG-alltype-firstGene_locustag.tsv
Rscript --vanilla $bin'R/TF-TG-alltype-firstGene_locustag_generator.R' -i $docs'TF-TG-alltype-firstGene.tsv' -a $docs'gene-bnumber.tsv' -o $docs'TF-TG-alltype-firstGene_locustag.tsv'

# 10) Generating TF-TG-alltype-firstGene_locustag_cop.tsv #eco.ope must be from RegulonDB
python $bin'python/TF-TG-alltype-firstGene_locustag_cop_generator.py' -p "$docs/eco.ope" -i "$docs/TF-TG-alltype-firstGene_locustag.tsv" -o "$docs/TF-TG-alltype-firstGene_locustag_cop.tsv"

# 11) Computing TP, TN, FN, FP and TPR and FPR according to the following definitions:

# TPR/specificity = TP / (TP + FN)
# FPR/sensivity = FP / (FP + TN)
# TP: Annotated binding site with a p-value < threshold
# FN: Annotated binding site with a p-value > threshold
# FP: Non-binding site with a p-value < threshold
# TN: Non-binding site with a p-value > threshold

cd $matrix_scan_regulonDB
for result in $(ls "$matrix_scan_regulonDB")
do
    TF=$(echo "$result" | cut -f1 -d "-")
    echo "$TF"
    grep -v ";" "$result" | cut -f1,2 | sed 's/|/\t/g' | sort | uniq | grep -v "#" > $TF'-predicted.tmp'
    Rscript --vanilla "$bin/R/Confusion_matrix_PSSM.R" -a $docs'TF-TG-alltype-firstGene_locustag_cop.tsv' -p $matrix_scan_regulonDB$TF'-predicted.tmp' -o $output_PSSM_RegulonDB$TF'-metrics.tsv' -t $TF -m "$output_PSSM_RegulonDB/FPR_TPR_table.tsv" -f "FALSE" && rm $TF'-predicted.tmp' 
done
