#Name:
#	merge_coverage_BLAST.R
#Author:
#	Hernandez Montserrat

#Version
#	v1
#Description
#	The script allows the merging of functional annotation
#	Input parameters:
#		--table1[directory]
#		--table2[directory]
#		--xMerge[character]
#		--yMerge[character]
#		--outdir[directory]
#
# 	Output
# 	1) A tsv file containing the merging annotation
# 
# Example:
# Rscript --vanilla merge_coverage_BLAST.R
# /home/emhernan/1_BBH_TFs/BLASTresults/ECaaq_RZaadb_blastP_b1_m8.tab
# /home/emhernan/1_BBH_TFs/BLASTresults/ECgene_len.tsv
# V1
# V1
# /home/emhernan/1_BBH_TFs/BLASTresults/ECaaq_RZaadb_blastN_b1_m8_v2.tab
# Rscript --vanilla merge_coverage_BLAST.R /home/emhernan/1_BBH_TFs/BLASTresults/ECaaq_RZaadb_blastP_b1_m8.tab /home/emhernan/1_BBH_TFs/BLASTresults/ECgene_len.tsv V1 V1 /home/emhernan/1_BBH_TFs/BLASTresults/ECaaq_RZaadb_blastN_b1_m8_v2.tab

# table1   <- "/home/emhernan/1_BBH_TFs/BLASTresults/ECaaq_RZaadb_blastP_b1_m8.tab"
# table2   <- "/home/emhernan/1_BBH_TFs/BLASTresults/ECgene_len.tsv"
# xMerge   <- "V1"
# yMerge   <- "V1"
# /home/emhernan/1_BBH_TFs/BLASTresults/ECaaq_RZaadb_blastN_b1_m8_v2.tab

# -*- encoding: utf-8 -*-

########################################## Merge coverage BLAST ########################################## 

print("...................................................Loading libraries!...................................................")

library(dplyr)

arg = commandArgs(trailingOnly = T)

if (length(arg)==0) {
  stop("Must supply input tables directories", call.=FALSE)
} else if (length(arg) == 1){
  stop("Enter input table 2 directory")
} else if (length(arg) == 2){
  stop("Enter the x column name to merge")
} else if (length(arg)==3){
  stop("Enter the y column name to merge")
} else if (length(arg)==4){
  stop("Enter the output directory")
}

table1   <- arg[1]
table2   <- arg[2]
xMerge   <- arg[3]
yMerge   <- arg[4]
outdir   <- arg[5]


print("...................................................Here are your input parameters!...................................................")
print(table1)
print(table2)
print(xMerge)
print(yMerge)
print(outdir)

print("...................................................Reading data...................................................")
p1 <- read.table(file = table1, header = FALSE, sep = "\t")
p2 <- read.table(file = table2, header = FALSE, sep = "\t")
print("...................................................Done!..................................................")

print("...................................................Processing data...................................................")

df <- left_join(p1, p2, by = c("V1" =  "V1")) %>%
		rename("qName" =  V1, "sName" =  V2.x, "peri" = V3, "alilen" = V4,"numMM" = V5,"nnGP" = V6,
			    "qSS" =  V7, "qSE" = V8,"sSS" = V9,"sSE" = V10,"Evalue" = V11,"bitScore" = V12,"qlen" = V2.y) %>%
		mutate(coveragePercent = ((qSE - qSS +1)*100)/qlen)


write.table(x = df, file = outdir, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
print("...................................................Done!..................................................")
