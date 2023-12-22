# Name:
#   section3-MainTable1.R
# Author:
#   Hernandez Benitez Ericka Montserrat
# Version:
#   v1
# Description:
#   The script computes the main table 1 with the following fields: TF_name, TF_locusTag, DNABindingRegion, ConservedRegion, effector_name, # of TGs
# Input parameters:
#   --networkPath[file path]
#   --TFsInfoPath[file path]
#   --motifEcoliPath[file path]
#   --motifRphaseoliPath[file path]
#   --outpath[file path]

# Output
#     1) A tsv file with the TFs' information from the TRN. 

# Rscript --vanilla section3-MainTable1.R 
# /space24/PGC/emhernan/3_TRN/output/oTF-TG-TRN-10-5-bgM1.tsv
# /space24/PGC/emhernan/3_TRN/OrthologousTFInfo/TF_Orthologous_table.tsv
# /space24/PGC/emhernan/3_TRN/OrthologousTFInfo/ECaaq_oTFs_motifs_info.tsv
# /space24/PGC/emhernan/3_TRN/OrthologousTFInfo/RZaaq_oTFs_motifs_info.tsv
# /space24/PGC/emhernan/3_TRN/tables/MainTab1.tsv


# Rscript --vanilla section3-MainTable1.R  /space24/PGC/emhernan/3_TRN/output/oTF-TG-TRN-10-5-bgM1.tsv /space24/PGC/emhernan/3_TRN/OrthologousTFInfo/TF_Orthologous_table.tsv /space24/PGC/emhernan/3_TRN/OrthologousTFInfo/ECaaq_oTFs_motifs_info.tsv /space24/PGC/emhernan/3_TRN/OrthologousTFInfo/RZaaq_oTFs_motifs_info.tsv /space24/PGC/emhernan/3_TRN/tables/MainTab1.tsv

# networkPath     <- "/space24/PGC/emhernan/3_TRN/output/oTF-TG-TRN-10-5-bgM1.tsv"
# TFsInfoPath     <- "/space24/PGC/emhernan/3_TRN/OrthologousTFInfo/TF_Orthologous_table.tsv"
# motifEcoliPath <- "/space24/PGC/emhernan/3_TRN/OrthologousTFInfo/ECaaq_oTFs_motifs_info.tsv"
# motifRphaseoliPath    <- "/space24/PGC/emhernan/3_TRN/OrthologousTFInfo/RZaaq_oTFs_motifs_info.tsv"
# outpath     <- "/space24/PGC/emhernan/3_TRN/tables/MainTab1.tsv"

############################ Loading libraries ############################
print("........................Loading libraries........................") 

library(dplyr)


arg = commandArgs(trailingOnly = T)

if (length(arg)==0) {
  stop("Must supply networkPath and TFsInfoPath file paths", call.=FALSE)
} else if (length(arg) == 1){
  stop("Enter TFsInfoPath file path")
} else if (length(arg) == 2){
  stop("Enter motifEcoliPath file path")
} else if (length(arg) == 3){
  stop("Enter motifRphaseoliPath file path")
} else if (length(arg) == 4){
  stop("Enter outpath file path")
} 


networkPath         <- arg[1]
TFsInfoPath         <- arg[2]
motifEcoliPath      <- arg[3]
motifRphaseoliPath  <- arg[4]
outpath             <- arg[5]


print("........................Here is your input data........................")

print(networkPath)
print(TFsInfoPath)
print(motifEcoliPath)
print(motifRphaseoliPath)
print(outpath)


print("........................Reading data........................")

network         <- read.table(file = networkPath, header =  TRUE, sep = "\t")
TFsinfo         <- read.table(file = TFsInfoPath, header =  TRUE, sep = "\t")
motifEcoli   <- read.table(file = motifEcoliPath, header =  TRUE, sep = "\t")
motifRphaseoli   <- read.table(file = motifRphaseoliPath, header =  TRUE, sep = "\t")

print("........................Processing data ........................")

df1 <- motifEcoli %>%
  filter(motifDesc == "DNA-Binding-Region" | motifDesc == "Conserved-Region") %>%
  filter(EC_locusTag %in% unique(network$TF_locusTag)) %>%
  rename(EcolimotifCoverage = motifCoverage) %>%
  rename(EcoliidentPercent = identPercent) %>%
  rename(mSSmotifEcoli = mSS) %>%
  rename(mSEmotifEcoli = mSE)

df2 <- motifRphaseoli  %>%
  filter(motifDesc == "DNA-Binding-Region" | motifDesc == "Conserved-Region") %>%
  filter(EC_locusTag %in% unique(network$TF_locusTag)) %>%
  rename(RphaseolimotifCoverage = motifCoverage) %>%
  rename(RphaseoliidentPercent = identPercent) %>%
  select(RphaseolimotifCoverage, RphaseoliidentPercent)

DNARegion <- bind_cols(df1, df2) %>%
                filter(motifDesc == "DNA-Binding-Region") %>%
                mutate(DNABindingRegion = paste(mSSmotifEcoli, mSEmotifEcoli, sep = "-")) %>%
                filter(RphaseolimotifCoverage == 100 & EcolimotifCoverage == 100 &
                       EcoliidentPercent >= 30 & RphaseoliidentPercent >= 30) %>%
            select(EC_locusTag, DNABindingRegion) %>%
            group_by(EC_locusTag) %>%
            summarise(DNABindingRegion = toString(DNABindingRegion)) %>%
            ungroup()

ConservedRegion <- bind_cols(df1, df2) %>%
                filter(motifDesc == "Conserved-Region") %>%
                mutate(ConservedRegion = paste(mSSmotifEcoli, mSEmotifEcoli, sep = "-")) %>%
                mutate(ConservedRegion_v = if_else(EcolimotifCoverage == 100 & 
                                  RphaseolimotifCoverage == 100 & 
                                  EcoliidentPercent >= 30 & 
                                 RphaseoliidentPercent >= 30, ConservedRegion, "NA")) %>%
                select(EC_locusTag, ConservedRegion_v)  %>%
                group_by(EC_locusTag) %>%
                summarise(ConservedRegion_v = toString(ConservedRegion_v)) %>%
                ungroup() %>%
                rename(ConservedRegion = ConservedRegion_v)
 
 

                
df3 <- network %>% 
  select(TF_name, TF_locusTag, TG_RZlocusTag) %>% 
  count(TF_name,TF_locusTag) %>%
  left_join((TFsinfo %>% select(EC_locusTag, effector_name) %>% distinct()), by = c("TF_locusTag" = "EC_locusTag")) %>%
  left_join(DNARegion, by = c("TF_locusTag" = "EC_locusTag")) %>%
  left_join(ConservedRegion, by = c("TF_locusTag" = "EC_locusTag")) %>%
  select(TF_name, TF_locusTag, DNABindingRegion, ConservedRegion, effector_name, n)  %>%
  rename("Number of TGs" = n) %>%
  rename("TF name" = TF_name) %>%
  rename("TF locus-tag" = TF_locusTag) %>%
  rename("Effector name" = effector_name) %>%
  rename("Conserved region" = ConservedRegion) %>%
  rename("DNA-binding region" = DNABindingRegion)

write.table(x = df3, file = outpath, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)