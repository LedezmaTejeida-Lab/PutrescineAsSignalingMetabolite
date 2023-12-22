#Name:
# oTF-TG-NetworkGenerator.R
#Author:
# Hernandez Montserrat
#Version
# v1
#Description
# The script allows the creation of a TF-TG tsv file 
# Input parameters:
#   --tgMatrixMotifPath[file path]
#   --MatrixMotifTfNamePath[file path]
#   --tfAnnotationPath[file path]
#   --tgTFRelationPath[file path]
#   --outpath[file path]
#
#   Output
#   1) A tsv file containt oTF-TG relations
# 
# Example:
# Rscript --vanilla oTF-TG-NetworkGenerator.R
# /space24/PGC/emhernan/3_TRN/output/p1_tmp
# /space24/PGC/emhernan/3_TRN/output/p2_tmp
# /space24/PGC/emhernan/3_TRN/output/p3_tmp
# /space24/PGC/emhernan/3_TRN/output/p4_tmp
# /space24/PGC/emhernan/3_TRN/output/p5_tmp
# /space24/PGC/emhernan/3_TRN/output/p5_tmp
# /space24/PGC/emhernan/3_TRN/output/TFOrthologous-TGInteractions.tsv
# Rscript --vanilla oTF-TG-NetworkGenerator.R /space24/PGC/emhernan/3_TRN/output/p1_tmp /space24/PGC/emhernan/3_TRN/output/p2_tmp /space24/PGC/emhernan/3_TRN/output/p3_tmp /space24/PGC/emhernan/3_TRN/output/p4_tmp /space24/PGC/emhernan/3_TRN/output/p5_tmp /space24/PGC/emhernan/3_TRN/output/p6_tmp /space24/PGC/emhernan/3_TRN/output/TFOrthologous-TGInteractions.tsv
# -*- encoding: utf-8 -*-

# tgMatrixMotifPath			<- "/space24/PGC/emhernan/3_TRN/output/p1_tmp"
# MatrixMotifTfNamePath		<- "/space24/PGC/emhernan/3_TRN/output/p2_tmp"
# tfAnnotationPath			<- "/space24/PGC/emhernan/3_TRN/output/p3_tmp"
# tgOrthologousAnnotationPath	<- "/space24/PGC/emhernan/3_TRN/output/p4_tmp"
# tgTFRelationPath			<- "/space24/PGC/emhernan/3_TRN/output/p5_tmp"
# tgRZlocustagRelationPath	<- "/space24/PGC/emhernan/3_TRN/output/p6_tmp"
# outpath						<- "/space24/PGC/emhernan/3_TRN/output/TFOrthologous-TGInteractions.tsv"

############################## Main program ###############################

############################ Loading libraries ############################

library(tidyverse)
library(fuzzyjoin)

arg = commandArgs(trailingOnly = T)

if (length(arg)==0) {
  stop("Must supply tgMatrixMotif and MatrixMotifTfName file paths", call.=FALSE)
} else if (length(arg) == 1){
  stop("Enter MatrixMotifTfName file path")
} else if (length(arg) == 2){
  stop("Enter tfAnnotation file path")
} else if (length(arg) == 3){
  stop("Enter tgOrthologousAnnotation file path")
} else if (length(arg) == 4){
  stop("Enter tgTFRelation file path")
} else if (length(arg) == 5){
  stop("Enter tgRZlocustagRelation file path")
} else if (length(arg) == 6){
  stop("Enter output file path")
}

tgMatrixMotifPath  			<- arg[1]
MatrixMotifTfNamePath 		<- arg[2]
tfAnnotationPath 			<- arg[3]
tgOrthologousAnnotationPath <- arg[4]
tgTFRelationPath 			<- arg[5]
tgRZlocustagRelationPath	<- arg[6]
outpath 					<- arg[7]

###########################################################################
############################## Inputs #####################################
###########################################################################

print("..........................................Here are your parameters:..........................................")
print(tgMatrixMotifPath)
print(MatrixMotifTfNamePath)
print(tfAnnotationPath)
print(tgOrthologousAnnotationPath)
print(tgTFRelationPath)
print(tgRZlocustagRelationPath)
print(outpath)
print("..........................................Done!..........................................")

###########################################################################
############################## Reading data ###############################
###########################################################################

print("..........................................Reading data..........................................")


tgMatrixMotif <- read.table(file=tgMatrixMotifPath, sep = "\t", header= FALSE)
MatrixMotifTfName <- read.table(file=MatrixMotifTfNamePath, sep = "\t", header= FALSE)
tfAnnotation <- read.table(file=tfAnnotationPath, sep = "\t", header= FALSE)
tgOrthologousAnnotation <- read.table(file=tgOrthologousAnnotationPath, sep = "\t", header= FALSE)
tgTFRelation <- read.table(file=tgTFRelationPath, sep = "\t", header= FALSE)
tgRZlocustagRelation <- read.table(file=tgRZlocustagRelationPath, sep = "\t", header= TRUE)
print("..........................................Done!..........................................")

###########################################################################
############################ Processing  data #############################
###########################################################################

print("..........................................Processing data..........................................")

df <- tgMatrixMotif %>%
	rename( TG_oldRZlocusTag = V1, TG_RZname = V2, scanType= V3, consensusMotif = V4) %>%
	left_join( MatrixMotifTfName %>% rename(consensusMotif = V1, TF_name = V2), by = "consensusMotif") %>%
	regex_join(tfAnnotation %>% rename( TF_locusTag = V1, TF_name_v = V2, TF_RdbID_v = V3, TF_RZlocusTag = V4, TF_oldRZlocusTag = V5), by = c("TF_name" = "TF_name_v"), mode = "left", ignore_case = TRUE) %>%
	left_join( tgOrthologousAnnotation %>% rename (TG_locusTag = V1, TG_name_v = V2, TG_RdbID_v= V3, TG_oldRZlocusTag_v = V4), by = c("TG_oldRZlocusTag" = "TG_oldRZlocusTag_v")) %>%
	tidyr::unite("key",TF_name,TG_RdbID_v, remove = FALSE, sep= "-", na.rm = TRUE) %>%
	left_join( tgTFRelation %>% rename (TF_RdbID = V1, TG_RdbID = V2, TG_name = V3, key = V4) %>% mutate(OrthologousInteration = TRUE), by = "key") %>%
	mutate(OrthologousInteration = ifelse(is.na(OrthologousInteration),FALSE, TRUE)) %>%
	left_join(tgRZlocustagRelation)	%>%
	select(scanType, consensusMotif, TF_RdbID_v, TF_name, TF_name_v, TF_locusTag, TF_RZlocusTag, TF_oldRZlocusTag, TG_RdbID_v, TG_name_v, TG_locusTag, TG_RZlocusTag, TG_oldRZlocusTag, TG_RZname, OrthologousInteration)
print("..........................................Done!..........................................")

###########################################################################
############################## Saving  data ###############################
###########################################################################

print("..........................................Saving data..........................................")
write.table(x = df, file = outpath ,sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
print("..........................................Done!..........................................")
