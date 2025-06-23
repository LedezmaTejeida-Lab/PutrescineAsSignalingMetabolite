# Name:
#   motifConservationAnalysis.R
# Author:
#   Hernandez Benitez Ericka Montserrat
# Version:
#   v1
# Description:
# The script computes the conservation of motifs 1 = conserved, 0 = no conserved, NA = no reported.
# The conservation analysis cosiders an deintity >= 40 and a coverage >=80 in both directions
# Input parameters:
#   --inpath[directory]
#   Output
#     1) A tsv file with the following columns: EC_locusTag, motifDesc, Motif_Conservation
#
# Rscript --vanilla motifConservationAnalysis.R
# /home/emhernan/2_MotifConservation/motifsInfo/
# Rscript --vanilla motifConservationAnalysis.R /home/emhernan/2_MotifConservation/motifsInfo/

# inpath <- "/home/emhernan/2_MotifConservation/motifsInfo/"


# -*- encoding: utf-8 -*-


arg = commandArgs(trailingOnly = T)
if (length(arg)==0) {
  stop("Must supply input directory", call.=FALSE)
}


inpath  <- arg[1]


print(".................................Loading libraries.................................")
library(dplyr)

print(".......................Reading data.......................")

Ecoli     <- read.table(file = paste(inpath, "ECaaq_oTFs_motifs_info.tsv", sep = ""), header = TRUE, sep = "\t")
Rphaseoli <- read.table(file = paste(inpath, "RZaaq_oTFs_motifs_info.tsv", sep = ""), header = TRUE, sep = "\t" )

print(".......................Processing data.......................")

df1 <- Ecoli %>%
  rename("EcolimotifCoverage" = motifCoverage) %>%
  rename("EcoliidentPercent" = identPercent) %>%
  na.omit()

df2 <- Rphaseoli %>%
  rename("RphaseolimotifCoverage" = motifCoverage) %>%
  rename("RphaseoliidentPercent" = identPercent) %>%
  select(RphaseolimotifCoverage, RphaseoliidentPercent) %>%
  na.omit()


df3 <-  Ecoli %>%
  select(EC_locusTag, motifDesc) %>%
  filter(is.na(motifDesc)) %>%
  mutate(motifDesc = NA) %>%
  mutate(Motif_Conservation = NA)

TFs <- bind_cols(df1, df2) %>%
  mutate(Motif_Conservation = if_else(EcoliidentPercent >= 40 & RphaseoliidentPercent >= 40 & EcolimotifCoverage >=80 & RphaseolimotifCoverage >= 80, 1, 0)) %>%
  mutate(motifDesc = sub("\\s+$", "", motifDesc)) %>%
  select(EC_locusTag, motifDesc, Motif_Conservation) %>%
  bind_rows(df3)

   



write.table(x = TFs, file = paste(inpath, "motifsConservationIdent40Coverage80.tsv", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
