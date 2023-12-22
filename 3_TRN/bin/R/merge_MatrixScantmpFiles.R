#Name:
#	merge_MatrixScantmpFiles.R
#Author:
#	Hernandez Montserrat
#Version
#	v1
#Description
#	The script allows the merging of columns from the RSAT program output of scan-matrix 
#	Input parameters:	     
#		--table1[file path]
#		--table2[file path]
#		--outdir[file path]
#
# 	Output
# 	1) A tsv file containing the merging tables
# 
# Example:
# Rscript --vanilla merge_MatrixScantmpFiles.R
# /space24/PGC/emhernan/3_TRN/output/scan-result-10-5-50-bgM1.tmp
# /space24/PGC/emhernan/3_TRN/output/relation_tmp
# /space24/PGC/emhernan/3_TRN/output/scan-result-10-5-50-bgM1.tmp2
# Rscript --vanilla merge_MatrixScantmpFiles.R /space24/PGC/emhernan/3_TRN/output/scan-result-10-5-50-bgM1.tmp /space24/PGC/emhernan/3_TRN/output/relation_tmp /space24/PGC/emhernan/3_TRN/output/scan-result-10-5-50-bgM1.tmp2


# -*- encoding: utf-8 -*-

# table1 <- "/space24/PGC/emhernan/3_TRN/output/scan-result-10-5-50-bgM1.tmp"
# table2 <- "/space24/PGC/emhernan/3_TRN/output/relation_tmp"
# outdir <- "/space24/PGC/emhernan/3_TRN/output/scan-result-10-5-50-bgM1.tmp2"

########################################## Merge Functional annotation ########################################## 


arg = commandArgs(trailingOnly = T)

if (length(arg)==0) {
  stop("Must supply table file paths", call.=FALSE)
} else if (length(arg) == 1){
  stop("Enter input table 2 file path")
} else if (length(arg)==2){
  stop("Enter the output file path")
} 

table1   <- arg[1]
table2   <- arg[2]
outdir   <- arg[3]

print("...................................................Here are your input parameters!...................................................")
library(tidyverse)

print(table1)
print(table2)
print(outdir)

print("...................................................Reading data...................................................")

p1 <- read.table(file = table1, header = FALSE, sep = "\t")
p2 <- read.table(file = table2, header = FALSE, sep = "\t")

print("...................................................Processing data...................................................")

df <- dplyr::left_join(p1, p2, by = c("V3" = "V1" )) %>%
		dplyr::select(V1, V2.x,V2.y,V4,V5,V6)
print("...................................................Done!..................................................")

print("...................................................Saving data...................................................")
write.table(x = df, file = outdir, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
print("...................................................Done!..................................................")

