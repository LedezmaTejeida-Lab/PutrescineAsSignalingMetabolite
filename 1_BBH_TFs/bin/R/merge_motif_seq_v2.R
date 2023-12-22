#	Name:
#		merge_motif_seq_v2.R
#	Author:
#		Hernandez Montserrat
#	Version
#		v2
#	Description
#		The script allows the merging of functional annotation
#	Input parameters:
#		--indir[directory]          
#		--outdir[output file path]
#		--gene_aa_seq[file name]     
#		--motifs_seq_relation[file name]   
# 		--prot_lt [file name]
#
# 	Output
# 		1) A tsv file containing the merging annotation
# 
#	Example:
#		Rscript --vanilla merge_merge_motif_seq_v2.R
#		/home/emhernan/1_BBH_TFs/motifsInfo/
#		/home/emhernan/1_BBH_TFs/motifsInfo/motifs_seq_relation_v3.tsv
#		gene_aa_seq.tsv
#		motifs_seq_relation_v2.tsv
#		prot_lt_tmp
#		Rscript --vanilla merge_motif_seq_v2.R /home/emhernan/1_BBH_TFs/motifsInfo/ /home/emhernan/1_BBH_TFs/motifsInfo/motifs_seq_relation_v3.tsv gene_aa_seq.tsv motifs_seq_relation_v2.tsv prot_lt_tmp
# -*- encoding: utf-8 -*-

# indir <- "/home/emhernan/1_BBH_TFs/motifsInfo/"
# outdir <- "/home/emhernan/1_BBH_TFs/motifsInfo/motifs_seq_relation_v3.tsv"
# gene_aa_seq <- "gene_aa_seq.tsv"
# motifs_seq_relation <- "motifs_seq_relation_v2.tsv"
# prot_lt <- "prot_lt_tmp"

########################################## Merge Functional annotation ########################################## 


arg = commandArgs(trailingOnly = T)

if (length(arg)==0) {
  stop("Must supply input and Output directories", call.=FALSE)
} else if (length(arg) == 1){
  stop("Enter output directory")
} else if (length(arg) == 2){
  stop("Enter aminoacid sequences file")
} else if (length(arg)== 3){
  stop("Enter motif sequence relation file")
} else if (length(arg)== 4){
  stop("Enter protein_id locustag relation file")
} 


indir				<- arg[1]
outdir				<- arg[2]
gene_aa_seq 		<- arg[3]
motifs_seq_relation	<- arg[4]
prot_lt				<- arg[5]





print("...................................................Here are your input parameters!...................................................")
print(indir)
print(outdir)
print(gene_aa_seq)
print(motifs_seq_relation)
print(prot_lt)



print("...................................................Reading data...................................................")
library(tidyr)
library(stringr)
library(dplyr)

setwd(indir)

p1 <- read.table(file = motifs_seq_relation, header = FALSE, sep = "\t")
p2 <- read.table(file = gene_aa_seq, header = FALSE, sep = "\t")
p3 <- read.table(file = prot_lt, header = TRUE, sep = "\t")


print("...................................................Done!..................................................")

print("...................................................Processing data...................................................")
p1 <- p1 %>% rename ("geneName" = V1, "locusTag" = V2, "motifDesc" = V3, "start" = V4, "end" = V5)

				
df <- p2 %>%
	rename("proteinID" = V1, "proteinSeq" = V2) %>%
	full_join(p3, by = c("proteinID"= "ProteinID")) %>%
	rename("locusTag" = Locus_tag) %>%
	right_join(p1, by = c( "locusTag" = "locusTag")) %>%
	rename("mSS" = start, "mSE" = end) %>%
	mutate(motifsSeq = str_sub(proteinSeq, mSS, mSE)) %>%
	select(locusTag, motifDesc, mSS, mSE, motifsSeq) %>%
	filter(!is.na(motifsSeq))


print("...................................................Done!..................................................")

print("...................................................Saving data...................................................")

write.table(x = df, file = outdir, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
print("...................................................Done!..................................................")