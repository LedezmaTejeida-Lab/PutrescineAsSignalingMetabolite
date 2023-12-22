#Name:
#	merge_Etables_p1_v3.R
#Author:
#	Hernandez Montserrat
#Version
#	v3
#Description
#	The script allows the merging of several df
#	Input parameters:
#		--indir[directory]
#		--outdir[file path]
#		--EC_BLASTid[file name]
#		--RZ_BLASTid[file name]
#		--EC_fant[file name]
#		--RZ_fant[file name]
#		--aaseq_EC[file path]
#		--BLAST_ECq[file path]
#
# 	Output
# 	1) A tsv file containing the merging annotation
# 
# Example:
# Rscript --vanilla merge_Etables_p1_v3.R
# /home/emhernan/1_BBH_TFs/tables/
# /home/emhernan/1_BBH_TFs/tables/Orthologous_table_p1.tsv
# Orthologous_ECaaq_RZaad_BLASTid_protein.tsv
# merge_annotation.tsv
# Orthologous_RZaad_ECaaq_BLASTid_protein.tsv
# Orthologous_RZ_annotation.tsv
# /home/emhernan/1_BBH_TFs/motifsInfo/gene_aa_seq.tsv
# /home/emhernan/1_BBH_TFs/tables/BLAST_ECq_tmp

# -*- encoding: utf-8 -*-

# Rscript --vanilla merge_Etables_p1_v3.R /home/emhernan/1_BBH_TFs/tables/ /home/emhernan/1_BBH_TFs/tables/Orthologous_table_p1.tsv Orthologous_ECaaq_RZaad_BLASTid_protein.tsv merge_annotation.tsv Orthologous_RZaad_ECaaq_BLASTid_protein.tsv Orthologous_RZ_annotation.tsv /home/emhernan/1_BBH_TFs/motifsInfo/gene_aa_seq.tsv /home/emhernan/1_BBH_TFs/tables/BLAST_ECq_tmp 

# indir   	<-  "/home/emhernan/1_BBH_TFs/tables/"
# outdir   	<-  "/home/emhernan/1_BBH_TFs/tables/Orthologous_table_p1.tsv"
# EC_BLASTid  <-  "Orthologous_ECaaq_RZaad_BLASTid_protein.tsv"
# EC_fant   	<-  "/home/emhernan/1_BBH_TFs/annotation/merge_annotation.tsv"
# RZ_BLASTid  <-  "Orthologous_RZaad_ECaaq_BLASTid_protein.tsv"
# RZ_fant 	<-  "Orthologous_RZ_annotation.tsv"
# aaseq_EC 	<-	"/home/emhernan/1_BBH_TFs/motifsInfo/gene_aa_seq.tsv"
# BLAST_ECq 	<-	"BLAST_ECq_tmp"


########################################## Merge Functional annotation ########################################## 


arg = commandArgs(trailingOnly = T)

if (length(arg)==0) {
  stop("Must supply input output directories", call.=FALSE)
} else if (length(arg)==1){
  stop("Enter output directory")
} else if (length(arg)==2){
  stop("Enter EC_BLASTid file name")
} else if (length(arg)==3){
  stop("Enter EC_fant file name")
} else if (length(arg)==4){
  stop("Enter RZ_BLASTid file name")
} else if (length(arg)==5){
  stop("Enter RZ_fant file name")
} else if (length(arg)==6){
  stop("Enter aaseq_EC file path")
} else if (length(arg)==7){
  stop("Enter BLAST_ECq file path")
} 

indir   	<- arg[1]
outdir   	<- arg[2]
EC_BLASTid  <- arg[3]
EC_fant   	<- arg[4]
RZ_BLASTid  <- arg[5]
RZ_fant 	<- arg[6]
aaseq_EC <- arg[7]
BLAST_ECq <- arg[8]

print("...................................................Here are your input parameters!...................................................")


print(indir)
print(outdir)
print(EC_BLASTid)
print(EC_fant)
print(RZ_BLASTid)
print(RZ_fant)
print(aaseq_EC)
print(BLAST_ECq)



print("...................................................Reading data...................................................")
library(tidyr)
library(stringr)
library(dplyr)

setwd(indir)
p1 <- as_tibble(read.table(file = EC_BLASTid, header = FALSE, sep = "\t"))
p2 <- as_tibble(read.table(file = EC_fant, header = TRUE, sep = "\t"))
p3 <- as_tibble(read.table(file = RZ_BLASTid, header = FALSE, sep = "\t"))
p4 <- as_tibble(read.table(file = RZ_fant, header = FALSE, sep = "\t"))
p5 <- as_tibble(read.table(file = aaseq_EC , header = FALSE, sep = "\t"))
p6 <- as_tibble(read.table(file = BLAST_ECq , header = TRUE, sep = "\t"))

print("...................................................Done!..................................................")

print("...................................................Processing data...................................................")

colnames(p1) <- c("ECBLAST_ID","RZBLAST_ID_inc","EC_proteinID")
colnames(p3) <- c("RZBLAST_ID","RZ_proteinID","RZBLAST_ID_inc_v")
colnames(p4) <- c("RZ_locusTag","RZ_Product","RZ_proteinID_v" )
colnames(p5) <- c("proteinID_v","Prot_seq")


df <- left_join(p1, p2, by = c("EC_proteinID" = "ProteinID")) %>%
		left_join(p3, by = c("RZBLAST_ID_inc" = "RZBLAST_ID_inc_v")) %>%
		left_join(p4, by = c("RZ_proteinID" = "RZ_proteinID_v")) %>%
		left_join(p5, by = c("EC_proteinID" = "proteinID_v")) %>%
		left_join(p6, by = c("ECBLAST_ID" = "qName")) %>%
		mutate(BLASTseq = str_sub(Prot_seq, qSS, qSE)) %>% 
		rename("EC_locusTag" = Locus_tag) %>%
		select(ECBLAST_ID, RZBLAST_ID, EC_locusTag, NCBI_name, Regulondb_name, 
			   Abasy_name, Ecocyc_name, Synonyms, RegulondbID, EC_proteinID, RZ_locusTag, 
			   RZ_proteinID, qSS, qSE, BLASTseq, EC_product, RZ_Product)


print("...................................................Done!..................................................")


print("...................................................Saving data...................................................")
write.table(x = df, file = outdir, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
print("...................................................Done!..................................................")
