#Name:
# oTF-TG-NetworkGenerator.R
#Author:
# Hernandez Montserrat
#Version
# v1
#Description
# The script allows the creation of a TF-TG tsv file 
# Input parameters:
#   --tf_tg_matrixscanPath[file path]
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
# /space24/PGC/emhernan/3_TRN/output/TFOrthologous-TGInteractions.tsv
# Rscript --vanilla oTF-TG-NetworkGenerator.R /space24/PGC/emhernan/3_TRN/output/p1_tmp /space24/PGC/emhernan/3_TRN/output/p2_tmp /space24/PGC/emhernan/3_TRN/output/p3_tmp /space24/PGC/emhernan/3_TRN/output/p4_tmp /space24/PGC/emhernan/3_TRN/output/p5_tmp /space24/PGC/emhernan/3_TRN/output/p6_tmp /space24/PGC/emhernan/3_TRN/output/TFOrthologous-TGInteractions.tsv
# -*- encoding: utf-8 -*-

# tf_tg_matrixscanPath				<- "/space24/PGC/emhernan/3_TRN/output/p1_tmp"
# tfAnnotationPath					<- "/space24/PGC/emhernan/3_TRN/output/p2_tmp"
# tgOrthologousAnnotationPath		<- "/space24/PGC/emhernan/3_TRN/output/p3_tmp"
# tgRZlocustagRelationPath			<- "/space24/PGC/emhernan/3_TRN/output/p4_tmp"
# tgTFRelationPath					<- "/space24/PGC/emhernan/3_TRN/output/p5_tmp"
# outpath							<- "/space24/PGC/emhernan/3_TRN/output/TFOrthologous-TGInteractions.tsv"

############################## Main program ###############################

############################ Loading libraries ############################

library(tidyverse)
library(fuzzyjoin)

arg = commandArgs(trailingOnly = T)

if (length(arg)==0) {
  stop("Must supply tf_tg_matrixscanPath and tfAnnotationPath file paths", call.=FALSE)
} else if (length(arg) == 1){
  stop("Enter tfAnnotationPath file path")
} else if (length(arg) == 2){
  stop("Enter tgOrthologousAnnotation file path")
} else if (length(arg) == 3){
  stop("Enter tgTFRelation file path")
} else if (length(arg) == 4){
  stop("Enter tgRZlocustagRelation file path")
} else if (length(arg) == 5){
  stop("Enter output file path")
}

tf_tg_matrixscanPath  			<- arg[1]
tfAnnotationPath 			<- arg[2]
tgOrthologousAnnotationPath <- arg[3]
tgRZlocustagRelationPath	<- arg[4]
tgTFRelationPath 			<- arg[5]
outpath 					<- arg[6]

###########################################################################
############################## Inputs #####################################
###########################################################################

print("..........................................Here are your parameters:..........................................")
print(tf_tg_matrixscanPath)
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


tf_tg_matrixscan <- read.table(file=tf_tg_matrixscanPath, sep = "\t", header= FALSE)
tfAnnotation <- read.table(file=tfAnnotationPath, sep = "\t", header= FALSE)
tgOrthologousAnnotation <- read.table(file=tgOrthologousAnnotationPath, sep = "\t", header= FALSE)
tgTFRelation <- read.table(file=tgTFRelationPath, sep = "\t", header= FALSE)
tgRZlocustagRelation <- read.table(file=tgRZlocustagRelationPath, sep = "\t", header= TRUE)
print("..........................................Done!..........................................")

###########################################################################
############################ Processing  data #############################
###########################################################################

print("..........................................Processing data..........................................")

df <- tf_tg_matrixscan %>%
	rename( TG_oldRZlocusTag = V1, TG_RZname = V2, scanType= V3, TFname = V4) %>%
	regex_join(tfAnnotation %>% rename( TF_locusTag = V1, TFname_v = V2, TF_RdbID = V3, TF_RZlocusTag = V4, TF_oldRZlocusTag = V5), by = c("TFname" = "TFname_v") , mode = "left", ignore_case = TRUE) %>%
	left_join(tgOrthologousAnnotation %>% rename (TG_locusTag = V1, TG_name = V2, TG_RdbID= V3, TG_oldRZlocusTag_v = V4), by = c("TG_oldRZlocusTag" = "TG_oldRZlocusTag_v")) %>%
	tidyr::unite("key",TFname,TG_RdbID, remove = FALSE, sep= "-", na.rm = TRUE) %>%
	left_join( tgTFRelation %>% rename (TF_RdbID_v = V1, TG_RdbID_v = V2, TG_name_v = V3, key = V4) %>% mutate(OrthologousInteration = TRUE), by = "key") %>%
	mutate(OrthologousInteration = ifelse(is.na(OrthologousInteration),FALSE, TRUE)) %>%
	left_join(tgRZlocustagRelation)	%>%
	select(scanType, TF_RdbID, TFname, TF_locusTag, TF_RZlocusTag, TF_oldRZlocusTag, TG_RdbID, TG_name, TG_locusTag, TG_RZlocusTag, TG_oldRZlocusTag, TG_RZname, OrthologousInteration)
print("..........................................Done!..........................................")

###########################################################################
############################## Saving  data ###############################
###########################################################################

print("..........................................Saving data..........................................")
write.table(x = df, file = outpath ,sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
print("..........................................Done!..........................................")
