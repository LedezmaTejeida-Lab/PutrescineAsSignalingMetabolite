#Name:
#	merge_Etables_p3_v5.R
#Author:
#	Hernandez Montserrat
#Version
#	v5
#Description
#	The script allows the merging of several df
#	Input parameters:
#		--Metabolites_relation[file path]
#		--RZ_relation_v2[file path]
#		--table_p2[file path]
#		--outdir[file path]

# 	Output
# 	1) A tsv file containing the merging annotation
# 	
# Example:
# Rscript --vanilla merge_Etables_p3_v5.R
# /home/emhernan/1_BBH_TFs/tables/Metabolites_relation_tmp
# /home/emhernan/1_BBH_TFs/annotation/RZ_functional_annotation_v2.tsv
# /home/emhernan/1_BBH_TFs/tables/Orthologous_table_p2.tsv
# /home/emhernan/1_BBH_TFs/tables/Orthologous_table.tsv

# Rscript --vanilla merge_Etables_p3_v5.R /home/emhernan/1_BBH_TFs/tables/Metabolites_relation_tmp /home/emhernan/1_BBH_TFs/annotation/RZ_functional_annotation_v2.tsv /home/emhernan/1_BBH_TFs/tables/Orthologous_table_p2.tsv /home/emhernan/1_BBH_TFs/tables/Orthologous_table.tsv


# -*- encoding: utf-8 -*-

# Metabolites_relation  <- "/home/emhernan/1_BBH_TFs/tables/Metabolites_relation_tmp"
# RZ_relation_v2      <- "/home/emhernan/1_BBH_TFs/annotation/RZ_functional_annotation_v2.tsv"
# table_p2        <- "/home/emhernan/1_BBH_TFs/tables/Orthologous_table_p2.tsv"
# outdir        <- "/home/emhernan/1_BBH_TFs/tables/Orthologous_table.tsv"

########################################## Merge Functional annotation ########################################## 


arg = commandArgs(trailingOnly = T)

if (length(arg)==0) {
  stop("Must supply Metabolites_relation RZ_relation_v2 path", call.=FALSE)
} else if (length(arg) == 1){
  stop("Enter RZ_relation_v2 path name")
} else if (length(arg) == 2){
  stop("Enter table_p2 path name")
} else if (length(arg) == 3){
  stop("Enter outdir path name")
}

Metabolites_relation    <- arg[1]
RZ_relation_v2          <- arg[2]
table_p2                <- arg[3]
outdir                  <- arg[4]
print("...................................................Here are your input parameters!...................................................")


print(Metabolites_relation)
print(RZ_relation_v2)
print(table_p2)
print(outdir)

print("...................................................Reading data...................................................")
library(tidyr)
library(stringr)
library(dplyr)
library(data.table)

p1 <- as_tibble(read.table(file = Metabolites_relation, header = TRUE, sep = "\t"))
p2 <- as_tibble(read.table(file = RZ_relation_v2, header = FALSE, sep = "\t"))
p3 <- as_tibble(read.table(file = table_p2, header = TRUE, sep = "\t"))

colnames(p2) <- c("RZ_locusTag_new", "RZ_locusTag_old")

print("...................................................Done!..................................................")

print("...................................................Processing data...................................................")

metabolites <- p1 %>%
              group_by(bnumber) %>%
              summarise(effector_name = toString(effector_name)) %>%
              ungroup() 


df <- left_join(p3, metabolites, by = c("EC_locusTag"="bnumber")) %>%
      left_join(p2, by = c("RZ_locusTag"="RZ_locusTag_new")) %>%
      select(ECBLAST_ID, RZBLAST_ID, EC_locusTag, NCBI_name, Regulondb_name, Abasy_name, Ecocyc_name, Synonyms, RegulondbID, EC_proteinID, RZ_locusTag, RZ_locusTag_old, RZ_proteinID, qSS, qSE, BLASTseq, EC_product, RZ_Product, motifDesc, mSS, mSE, motifsSeq, effector_name)



print("...................................................Saving data...................................................")
write.table(x = df, file = outdir, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
print("...................................................Done!..................................................")


