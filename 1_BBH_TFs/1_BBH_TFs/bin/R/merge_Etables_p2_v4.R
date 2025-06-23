#Name:
#	merge_Etables_p2_v4.R
#Author:
#	Hernandez Montserrat
#Version
#	v4
#Description
#	The script allows the merging of several df
#	Input parameters:
#		--motif_relation[file path TF/Orthologous]
#		--table_p1[file path]
#		--outdir[output path]
#
# 	Output
# 	1) A tsv file containing the merging annotation
# 
# Example:
# Rscript --vanilla merge_Etables_p2_v4.R
# /home/emhernan/1_BBH_TFs/motifs_seq_relation_v3.tsv
# /home/emhernan/1_BBH_TFs/tables/Orthologous_table_p1.tsv
# /home/emhernan/1_BBH_TFs/tables/Orthologous_table_p2.tsv
# Rscript --vanilla merge_Etables_p2_v3.R /home/emhernan/1_BBH_TFs/motifs_seq_relation_v3.tsv /home/emhernan/1_BBH_TFs/tables/Orthologous_table_p1.tsv /home/emhernan/1_BBH_TFs/tables/Orthologous_table_p2.tsv

# motif_relation 	<- "/home/emhernan/1_BBH_TFs/motifsInfo/motifs_seq_relation_v3.tsv"
# table_p1 		<- "/home/emhernan/1_BBH_TFs/tables/Orthologous_table_p1.tsv"
# outdir 			<- "/home/emhernan/1_BBH_TFs/tables/Orthologous_table_p2.tsv"

# -*- encoding: utf-8 -*-


########################################## Merge Functional annotation ########################################## 


arg = commandArgs(trailingOnly = T)

if (length(arg)==0) {
  stop("Must supply motif_relation and table_p1 path files", call.=FALSE)
} else if (length(arg) == 1){
  stop("Enter table_p1 path name")
} else if (length(arg) == 2){
  stop("Enter output path")
} 


motif_relation 	<- arg[1]
table_p1 		<- arg[2]
outdir 			<- arg[3]

print("...................................................Here are your input parameters!...................................................")

print(motif_relation)
print(table_p1)
print(outdir)

library(tidyr)
library(stringr)
library(dplyr)

print("...................................................Done!..................................................")


print("...................................................Reading data...................................................")
p1 <- as_tibble(read.table(file = motif_relation, header = FALSE, sep = "\t"))
p2 <- as_tibble(read.table(file = table_p1, header = TRUE, sep = "\t"))

print("...................................................Done!..................................................")

print("...................................................Processing data...................................................")

df <-  p1 %>%
		rename("locusTag" = V1, "motifDesc" = V2 , "mSS" = V3, "mSE" = V4, "motifsSeq" = V5) %>%
		right_join(p2, by = c("locusTag" = "EC_locusTag")) %>% 
		rename("EC_locusTag" = locusTag) %>%
		select(ECBLAST_ID, RZBLAST_ID, EC_locusTag, NCBI_name, Regulondb_name, Abasy_name, Ecocyc_name, Synonyms, RegulondbID, EC_proteinID, RZ_locusTag, RZ_proteinID, qSS, qSE, BLASTseq, EC_product, RZ_Product, motifDesc, mSS, mSE, motifsSeq)


print("...................................................Saving data...................................................")
write.table(x = df, file = outdir, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
print("...................................................Done!..................................................")

