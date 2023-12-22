#########################################################################################################
######################################### Working directory #############################################
#########################################################################################################

cd /home/emhernan/1_BBH_TFs/

#########################################################################################################
##################################### Saving all work directories #####################################
#########################################################################################################

indir=/home/emhernan/5_BBH_PuuR/
bin=/home/emhernan/5_BBH_PuuR/bin/
docs=/home/emhernan/5_BBH_PuuR/docs/
Rhizobia=/home/emhernan/5_BBH_PuuR/docs/Rhizobia/
RhizobiumPhaseoli=/home/emhernan/5_BBH_PuuR/docs/RhizobiumPhaseoli/
output=/home/emhernan/5_BBH_PuuR/output/
multiple_aligment=/home/emhernan/5_BBH_PuuR/docs/multiple_aligment/
png=/home/emhernan/5_BBH_PuuR/png/


cd $indir

#########################################################################################################
############################################### BLAST BBH ###############################################
#########################################################################################################

## 0) Creating folders

mkdir $docs
mkdir $Rhizobia
mkdir $RhizobiumPhaseoli
mkdir $output
mkdir $png



## 1) Unziping files
cd $Rhizobia

for file in $(ls "$Rhizobia")
do
    unzip "$file" -d "${file%.zip}" *faa *fna *gbff
    mv "${file%.zip}/ncbi_dataset/data/"* "${file%.zip}"
    rm -r "${file%.zip}/ncbi_dataset"
    rm -r "$file"

done

cd $RhizobiumPhaseoli

for file in $(ls "$RhizobiumPhaseoli")
do
    unzip "$file" -d "${file%.zip}" *faa *fna *gbff
    mv "${file%.zip}/ncbi_dataset/data/"* "${file%.zip}"
    rm -r "${file%.zip}/ncbi_dataset"
    rm -r "$file"

done


## 2) Generate a proper fasta header for db indexing

cd "$Rhizobia"

for dir in $(ls "$Rhizobia")
do
    GCF=$(ls "$dir" | grep GCF_*)
    cat "$Rhizobia$dir/$GCF/protein.faa" | perl -nae 'if(/>(.+\.\d)(.+)\s\[(.+)\]/) {$c++; print ">gnl|Rhadb|$c|$1$2\n"} else {print "$_"}' > "$Rhizobia$dir/$GCF/protein.faaed"
done

cd "$RhizobiumPhaseoli"

for dir in $(ls "$RhizobiumPhaseoli")
do
    GCF=$(ls "$dir" | grep GCF_*)
    cat "$RhizobiumPhaseoli$dir/$GCF/protein.faa" | perl -nae 'if(/>(.+\.\d)(.+)\s\[(.+)\]/) {$c++; print ">gnl|Rhpdb|$c|$1$2\n"} else {print "$_"}' > "$RhizobiumPhaseoli$dir/$GCF/protein.faaed"
done



## 3) Run formatdb, generating an indexed DB
cd $Rhizobia

for dir in $(ls "$Rhizobia")
do
    GCF=$(ls "$dir" | grep GCF_*)
    mkdir "$Rhizobia$dir/$GCF/BLAST"
    formatdb -i "$Rhizobia$dir/$GCF/protein.faaed" -p T -o T
    for file in $(ls "$Rhizobia$dir/$GCF/protein.faaed"*)
    do
        mv "$file" "$Rhizobia$dir/$GCF/BLAST"
    done
done


cd $RhizobiumPhaseoli

for dir in $(ls "$RhizobiumPhaseoli")
do
    GCF=$(ls "$dir" | grep GCF_*)
    mkdir "$RhizobiumPhaseoli$dir/$GCF/BLAST"
    formatdb -i "$RhizobiumPhaseoli$dir/$GCF/protein.faaed" -p T -o T
    for file in $(ls "$RhizobiumPhaseoli$dir/$GCF/protein.faaed"*)
    do
        mv "$file" "$RhizobiumPhaseoli$dir/$GCF/BLAST"
    done
done


## 4) Running BLAST

cd $Rhizobia

for dir in $(ls "$Rhizobia")
do
    GCF=$(ls "$dir" | grep GCF_*)
    for file in $(ls "$Rhizobia$dir/$GCF/BLAST/protein.faaed")
    do
        blastall -p blastp -i "$file" -d "$RhizobiumPhaseoli"*"/GCF"*"/BLAST/protein.faaed" -b 1 -a 6 -e 0.001 -m 8 > "$Rhizobia$dir/$GCF/BLAST/Rhaaq_Rpaadb_blastP_b1_m8.tab"
        blastall -p blastp -i "$RhizobiumPhaseoli"*"/GCF"*"/BLAST/protein.faaed" -d "$file" -b 1 -a 6 -e 0.001 -m 8 > "$Rhizobia$dir/$GCF/BLAST/Rpaaq_Rhaadb_blastP_b1_m8.tab"
    done
done



## 5) Searching fo orthologous
## gnl|RZaadb|5641|WP_029531806.1 --> PuuR gnl|ECaadb|1189|NP_415815.1 --> gnl|Rhpdb|5559|WP_029531806.1

cd $Rhizobia

for dir in $(ls "$Rhizobia")
do
    GCF=$(ls "$dir" | grep GCF_*)
    Rpdb=$(paste <(awk '{if ($11 <= 0.001 && $3 >=30 && $12 >=50) {print }}' "$Rhizobia$dir/$GCF/BLAST/Rpaaq_Rhaadb_blastP_b1_m8.tab"| grep -w "gnl|Rhpdb|5559|WP_029531806.1" | cut -f1 | cut -f1-3 -d "|") <(awk '{if ($11 <= 0.001 && $3 >=30 && $12 >=50) {print }}' "$Rhizobia$dir/$GCF/BLAST/Rpaaq_Rhaadb_blastP_b1_m8.tab"| grep -w "gnl|Rhpdb|5559|WP_029531806.1" | cut -f2))
    Rhdb=$(paste <(awk '{if ($11 <= 0.001 && $3 >=30 && $12 >=50) {print }}' "$Rhizobia$dir/$GCF/BLAST/Rhaaq_Rpaadb_blastP_b1_m8.tab" | grep -w "gnl|Rhpdb|5559" | cut -f1 | cut -f1-3 -d "|") <(awk '{if ($11 <= 0.001 && $3 >=30 && $12 >=50) {print }}' "$Rhizobia$dir/$GCF/BLAST/Rhaaq_Rpaadb_blastP_b1_m8.tab" | grep -w "gnl|Rhpdb|5559" |cut -f2))
    grep -wf <(echo "$Rpdb" | perl -nae 'print "$F[1]\t$F[0]\n"') <(echo "$Rhdb") > "$Rhizobia$dir/$GCF/BLAST/Orthologous_PuuRID.tsv"
    grep -wf <(cut -f1 "$Rhizobia$dir/$GCF/BLAST/Orthologous_PuuRID.tsv") "$Rhizobia$dir/$GCF/BLAST/Rpaaq_Rhaadb_blastP_b1_m8.tab" > "$Rhizobia$dir/$GCF/BLAST/Rpaaq_RhaadbResult.tsv"
	grep -f <(cut -f1 "$Rhizobia$dir/$GCF/BLAST/Orthologous_PuuRID.tsv") <(grep -wf <(cut -f2 "$Rhizobia$dir/$GCF/BLAST/Orthologous_PuuRID.tsv") "$Rhizobia$dir/$GCF/BLAST/Rhaaq_Rpaadb_blastP_b1_m8.tab") > "$Rhizobia$dir/$GCF/BLAST/Rhaaq_RpaadbResult.tsv"

done




## 6) Generating only one report per query

cd $Rhizobia
Rpaaq_RhaadbResult=""
Rhaaq_RpaadbResult=""

for dir in $(ls "$Rhizobia")
do
    GCF=$(ls "$dir" | grep GCF_*)
    for file in $(ls "$Rhizobia$dir/$GCF/BLAST/Rpaaq_RhaadbResult.tsv")
    do
        Rpaaq_RhaadbResult=$(cat <(echo "$Rpaaq_RhaadbResult") "$file" )
    done
done



for dir in $(ls "$Rhizobia")
do
    GCF=$(ls "$dir" | grep GCF_*)
    for file in $(ls "$Rhizobia$dir/$GCF/BLAST/Rhaaq_RpaadbResult.tsv")
    do
        Rhaaq_RpaadbResult=$(cat <(echo "$Rhaaq_RpaadbResult") "$file")
    done
done

mkdir $output

cd $output
echo "$Rpaaq_RhaadbResult" > $output'Rpaaq_Rhaadb_PuuROrthologous.tsv'
echo "$Rhaaq_RpaadbResult" > $output'Rhaaq_Rpaadb_PuuROrthologous.tsv'

sed -i '1d' $output'Rpaaq_Rhaadb_PuuROrthologous.tsv'
sed -i '1d' $output'Rhaaq_Rpaadb_PuuROrthologous.tsv'

# 7) Plotting Main Figure 6 and Supp Figure 5
Rscript --vanilla $bin'R/section5-MainFigure6SuppFig5.R' $output $png



#8) Performing multiple aligment with clustalw

cd $multiple_aligment
clustalw -infile=$multiple_aligment'PuuR-BBH-sequences.faa' -outfile=$output'PuuR-BBH-sequences_cluAln.nex' -output=nexus -pwgapopen=14 -pwgapext=0.2 -gapopen=14 -gapext=0.2

# 9) Performing MP

cd $output

paup

[weights 2:2stpos; ctype 2_1:all;]

Begin paup;
	set autoclose=yes warntree=no warnreset=no;
	log start file=PuuR-BBH.log replace;
	exe PuuR-BBH-sequences_cluAln.nex;
	set criterion=parsimony; 
	cstatus;
	hsearch nbest=5;
	outgroup Escherichia_coli_NP_415815.1;
	Showtree all;
	describetrees 4/ plot=phylogram brlens=yes rootmethod=outgroup;
	SAVETREES file=PuuR-BBH.tree brlens=yes format = newick replace;
	bootstrap treeFile=MP_BootstrapTreeFilea.tre nreps=1000 conLevel=50 brLens=yes replace=yes format=Newick search=heuristic keepAll=yes/nbest=0;
	savetrees from=1 to=1 file=bootstrap_trees.tre format=newick savebootp=nodelabels;

