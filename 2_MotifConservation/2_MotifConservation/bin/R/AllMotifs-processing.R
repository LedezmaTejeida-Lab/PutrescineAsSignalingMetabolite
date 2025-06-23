# Name:
#   AllMotifs-processing.R
# Author:
#   Hernandez Benitez Ericka Montserrat
# Version:
#   v1
# Description:
#   The script computes the identity and coverage of the motifs reported in Ecocyc for each TF. The identity and coverage
#   are based on the BLASTp aligments of all the protein sequence. 
# Input parameters:
#   --inpath[directory]
#   --outpath[directory]
#   --queryFile[file name]
#   --seqIdentFile_dash[file name]
#   --seqIdentFile_noDash [file name] 
#   --queryString [query identifier in the BLASTp]
#   --outFile [filene name]
# Output
#     1) A tsv file with the following columns: NCBI_name, EC_locusTag, motifCoverage, identPercent, motifDesc
#
# Rscript --vanilla AllMotifs-processing.R
# /home/emhernan/2_MotifConservation/formated/
# /home/emhernan/2_MotifConservation/motifsInfo/
# ECaaq-queID.tsv
# ECaaq_RZaadb_blastP_seq_formated.tsv
# ECaaq_RZaadb_blastP_seq_formated_v2.tsv
# ECBLAST_ID
# ECaaq_oTFs_motifs_info.tsv
# Rscript --vanilla AllMotifs-processing.R /home/emhernan/2_MotifConservation/formated/ /home/emhernan/2_MotifConservation/motifsInfo/ ECaaq-queID.tsv ECaaq_RZaadb_blastP_seq_formated.tsv ECaaq_RZaadb_blastP_seq_formated_v2.tsv ECBLAST_ID ECaaq_oTFs_motifs_info.tsv

# inpath  <- "/home/emhernan/2_MotifConservation/formated/"
# outpath <- "/home/emhernan/2_MotifConservation/motifsInfo/"
# queryFile     <- "ECaaq-queID.tsv"
# seqIdentFile_dash  <- "ECaaq_RZaadb_blastP_seq_formated.tsv"
# seqIdentFile_noDash  <- "ECaaq_RZaadb_blastP_seq_formated_v2.tsv"
# queryString   <- "ECBLAST_ID"
# outFile       <- "ECaaq_oTFs_motifs_info.tsv"
# -*- encoding: utf-8 -*-


print(".................................Loading libraries.................................")

library(dplyr)
library(stringr)
library(DescTools)
library(stringi)

arg = commandArgs(trailingOnly = T)
if (length(arg)==0) {
  stop("Must supply input and output directories", call.=FALSE)
} else if (length(arg) == 1){
  stop("Enter output directory")
} else if (length(arg) == 2){
  stop("Enter query IDs file name")
} else if (length(arg)==3){
  stop("Enter sequence identity (dash) file name")
} else if (length(arg)==4){
  stop("Enter sequence identity (non-dash) file name")
} else if (length(arg)==5){
  stop("Enter query string")
} else if (length(arg)==6){
  stop("Enter output file name")
}

inpath        <- arg[1]
outpath       <- arg[2]
queryFile     <- arg[3]
seqIdentFile_dash  <- arg[4]
seqIdentFile_noDash  <- arg[5]
queryString   <- arg[6]
outFile       <- arg[7]




print(".................................Reading data.................................")

qID       <- read.table(file = paste(inpath, queryFile, sep = ""), header = FALSE, sep = "\t")
seqIdent_dash  <- read.table(file = paste(inpath, seqIdentFile_dash, sep = ""), header = TRUE, sep = "\t")
seqIdent_noDash  <- read.table(file = paste(inpath, seqIdentFile_noDash, sep = ""), header = TRUE, sep = "\t")
TFsInfo   <- read.table(file = paste(inpath, "TF_Orthologous_table.tsv", sep = ""), header = TRUE, sep = "\t")


print(".................................Processing data.................................")


colnames(qID) <- queryString
colnames(seqIdent_dash) <- c("start", "end", "seq_dash", "seqIdent_dash")
colnames(seqIdent_noDash) <- c("start", "end", "seq_noDash", "seqIdent_noDash")

df1 <- TFsInfo %>%
  select(all_of(queryString), EC_locusTag, NCBI_name, mSS, mSE, motifDesc, motifsSeq) %>%
  distinct() %>%
  na.omit() %>%
  as_tibble()
  
df1_na <- TFsInfo %>%
  select(NCBI_name, EC_locusTag, mSS, mSE) %>%
  filter(is.na(mSS) & is.na(mSE)) %>%
  distinct() %>%
  as_tibble()

df1_na <- df1_na %>%
            mutate(motifCoverage = NA) %>%
            mutate(identPercent = NA) %>%
            mutate(motifDesc = NA) %>%
            select(NCBI_name, EC_locusTag, motifCoverage, identPercent, motifDesc)

df2 <- bind_cols(qID, seqIdent_noDash) %>% 
  bind_cols((seqIdent_dash %>% select(seq_dash, seqIdent_dash))) %>%
  as_tibble()

df3 <- df1 %>%
  left_join(df2) %>%
  mutate(motifLen = str_length(motifsSeq))


motif <- mapply(seq, df3$mSS, df3$mSE)
BLAST <- mapply(seq, df3$start, df3$end)
df3$overlap <- mapply(Overlap, x = motif, y = BLAST)

df3 <- df3 %>%
  mutate(overlap = ifelse(overlap > 0, overlap + 1 , 0)) %>%
  mutate(motifCoverage = round((overlap * 100)/motifLen,2)) %>%
  mutate(motif_noDash = substr(seq_noDash, (mSS-start+1), (mSE-start+1))) %>%
  select(EC_locusTag, NCBI_name, motifCoverage, seq_dash, motif_noDash, seqIdent_dash, mSS, mSE, motifDesc)

motif_noDash_v <- df3$motif_noDash
motif_collapsed <- mapply(str_split, motif_noDash_v, "", SIMPLIFY = TRUE, USE.NAMES = FALSE)
df3$motif_pattern <- mapply(paste, motif_collapsed, collapse = '-*')

df3 <- df3 %>%
  mutate(motif_pattern = if_else(motif_pattern == "", "1", motif_pattern)) %>%
  mutate(motif_dash = str_extract(string = seq_dash, pattern = motif_pattern)) %>%
  select(EC_locusTag, NCBI_name, motifCoverage, motif_dash, motif_pattern, seq_dash, seqIdent_dash, mSS, mSE, motifDesc)

seq_dash_v <- df3$seq_dash
motif_pattern_v <- df3$motif_pattern
coordinates <- as.data.frame(t(mapply(str_locate, seq_dash_v, motif_pattern_v, SIMPLIFY = TRUE, USE.NAMES = FALSE)))
colnames(coordinates) <- c("Start_motif_dash", "End_motif_dash")


df3 <- bind_cols(df3, coordinates) %>%
    mutate(identityMotif_dash = substr(seqIdent_dash, Start_motif_dash, End_motif_dash)) %>%
    mutate(identPercent = round (((stri_length(identityMotif_dash) - stri_count_fixed(identityMotif_dash, " ") - stri_count_fixed(identityMotif_dash, "+")) * 100)/stri_length(identityMotif_dash), 2)) %>%
    mutate(motifDesc = sub("\\s+$", "", motifDesc)) %>%
    select(NCBI_name, EC_locusTag, motifCoverage, identPercent, motifDesc, mSS, mSE) %>%
    replace(is.na(.), 0)
    


df4 <- bind_rows(df3, df1_na)



print(".................................Saving data.................................")

write.table(x = df4, file = paste(outpath, outFile, sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")

