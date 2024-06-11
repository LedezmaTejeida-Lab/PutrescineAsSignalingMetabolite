# Name:
#   DNAbindingFilterTFs100Cov30Ident.R
# Author:
#   Hernandez Benitez Ericka Montserrat
# Version:
#   v1
# Description:
# The script computes the filter of TFs that have a DNA-binding motif with coverage = 100 and identoty >= 30 in both directions
# Input parameters:
#   --inpath[directory]
#   --outfile[file name]
#   Output
#     1) A tsv file with the following columns: EC_locusTag, NCBI_name, Regulondb_name, Abasy_name, Ecocyc_name, Synonyms
#
# Rscript --vanilla DNAbindingFilterTFs100Cov30Ident.R
# /space24/PGC/emhernan/3_TRN/OrthologousTFInfo/
# TF_Orthologous_table_filter100Cov40Iden.tsv
# Rscript --vanilla DNAbindingFilterTFs100Cov30Ident.R /space24/PGC/emhernan/3_TRN/OrthologousTFInfo/ TF_Orthologous_table_filter100Cov40Iden.tsv

# inpath <- "/space24/PGC/emhernan/3_TRN/OrthologousTFInfo/"
# outfile <- "TF_Orthologous_table_filter100Cov40Iden.tsv"

# -*- encoding: utf-8 -*-


arg = commandArgs(trailingOnly = T)
if (length(arg)==0) {
  stop("Must supply input directory and output file", call.=FALSE)
} else if (length(arg)==1) {
  stop("Must supply output file", call.=FALSE)
}


inpath <- arg[1]
outfile <- arg[2]

print(".................................Loading libraries.................................")

library(dplyr)


print(".......................Reading data.......................")
ecoli <- read.table(file = paste(inpath, "ECaaq_oTFs_DNAbinding_motifs_info.tsv", sep = ""), header = FALSE, sep = "\t")
rphaseoli <- read.table(file = paste(inpath, "RZaaq_oTFs_DNAbinding_motifs_info.tsv", sep = ""), header = FALSE, sep = "\t")
TF_orthologous_table <- read.table(file = paste(inpath, "TF_Orthologous_table.tsv", sep = ""), header = TRUE, sep = "\t")



print(".......................Processing data.......................")
ecoli <- ecoli %>% 
				rename( "geneName" = V1, "locusTag" = V2, "coverageEcoli" = V3, "identityEcoli" = V4, "motifDesc" = V5) %>%
				select(geneName, locusTag, coverageEcoli, identityEcoli)


rphaseoli <- rphaseoli %>% 
					rename( "geneName" = V1, "locusTag" = V2, "coverageRphaseoli" = V3, "identityRphaseoli" = V4, "motifDesc" = V5) %>%
					select(coverageRphaseoli, identityRphaseoli)

TFs_total <- bind_cols(ecoli, rphaseoli) %>%
			select(locusTag) %>%
			pull() %>%
			unique() %>%
			length()

df <- bind_cols(ecoli, rphaseoli) %>%
	mutate(strict_cutOff = if_else(identityEcoli >= 40 & identityRphaseoli >= 40 & coverageEcoli == 100 & coverageRphaseoli == 100, 1 , 0))


TFs_1 <- df %>%
			filter(strict_cutOff == 1) %>%
			select(locusTag) %>%
			pull() %>%
			unique()

TFs_0 <- df %>%
			filter(strict_cutOff == 0) %>%
			select(locusTag) %>%
			pull() %>%
			unique()

# unique(setdiff(TFs_0, TFs_1))
# unique(intersect(TFs_1, TFs_0))


IDsWithCutOff <- unique(setdiff(TFs_1, TFs_0))

TF_orthologous_tableFilter <- TF_orthologous_table %>%
				filter(EC_locusTag %in% IDsWithCutOff) %>%
				select(EC_locusTag, NCBI_name, Regulondb_name, Abasy_name, Ecocyc_name, Synonyms) %>%
				distinct()

print("......................Saving data.......................")

write.table(x = TF_orthologous_tableFilter,  file= paste(inpath, outfile, sep = ""),
			sep = "\t", row.names = FALSE, quote = FALSE)