## Main results from section 2 of article

library(dplyr)

ecoli <- read.table(file = "/home/emhernan/2_MotifConservation/motifsInfo/ECaaq_oTFs_motifs_info.tsv", header = TRUE, sep = "\t")
rphaseoli <- read.table(file = "/home/emhernan/2_MotifConservation/motifsInfo/RZaaq_oTFs_motifs_info.tsv", header = TRUE, sep = "\t")
conservation <-  read.table(file = "/home/emhernan/2_MotifConservation/motifsInfo/motifsConservationIdent40Coverage80.tsv", header = TRUE, sep = "\t")

#### TFs with NO reported motifs #####

print("Number of TFs with NO reported motifs Ecoli")


TFsWithoutMotifEcoli <- ecoli %>%
						filter(is.na(motifDesc)) %>%
						select(NCBI_name) %>%
						pull()

print(length(TFsWithoutMotifEcoli))
print(TFsWithoutMotifEcoli)

print("Number of TFs with NO reported motifs Rphaseoli")


TFsWithoutMotifRphaseoli <- rphaseoli %>%
						filter(is.na(motifDesc)) %>%
						select(NCBI_name) %>%
						pull()

print(length(TFsWithoutMotifRphaseoli))
print(TFsWithoutMotifRphaseoli)


#### TFs with reported motifs #####

print("Number of TFs with reported motifs Ecoli")


TFsWithMotifEcoli <- ecoli %>%
						filter(!is.na(motifDesc)) %>%
						select(NCBI_name) %>%
						pull() %>%
						unique()
print(length(TFsWithMotifEcoli))
print(TFsWithMotifEcoli)

print("Number of TFs with reported motifs Rphaseoli")


TFsWithMotifRphaseoli <- rphaseoli %>%
						filter(!is.na(motifDesc)) %>%
						select(NCBI_name) %>%
						pull() %>%
						unique()

print(length(TFsWithMotifRphaseoli))
print(TFsWithMotifRphaseoli)


#### How many motifs are reported for the 49 TFs? #####

numberTFs <- as.character(length(TFsWithMotifRphaseoli))

print(paste("Number of motifs annotated for", numberTFs, "TFs Ecoli"), sep = "")


motifsEcoli <- ecoli %>% 
				filter(!is.na(motifDesc)) %>%
				select(motifDesc) %>%
				pull()

print(length(motifsEcoli))
print(motifsEcoli)


print(paste("Number of motifs annotated for", numberTFs, "TFs Rphaseoli"), sep = "")


motifsRphaseoli <- rphaseoli %>% 
				filter(!is.na(motifDesc)) %>%
				select(motifDesc) %>%
				pull()

print(length(motifsRphaseoli))
print(motifsRphaseoli)


#### How many motifs on average per TF are reported? #####

print("Number of reported motifs on avergae per TF Ecoli")


motifsPerTF <- ecoli %>% 
				filter(!is.na(motifDesc)) %>%
				count(NCBI_name)

print(mean(motifsPerTF$n))
print(motifsPerTF)


print("Number of reported motifs on avergae per TF Rphaseoli")


motifsPerTFRphaseoli <- rphaseoli %>% 
				filter(!is.na(motifDesc)) %>%
				count(NCBI_name)

print(mean(motifsPerTFRphaseoli$n))
print(motifsPerTFRphaseoli)


#### How many TFs have at least one DNA-binding motif reported? #####

print("Number of TFs with at least one DNA-binding motif reported TF Ecoli")


TFsWithDNAbindingMotif <- ecoli %>% 
				filter(motifDesc == "DNA-Binding-Region") %>%
				select(NCBI_name) %>%
				unique() %>%
				pull()


print(length(TFsWithDNAbindingMotif))
print(TFsWithDNAbindingMotif)




print("Number of TFs with at least one DNA-binding motif reported TF Rphaseoli")


TFsWithDNAbindingMotifRphaseoli <- rphaseoli %>% 
				filter(motifDesc == "DNA-Binding-Region") %>%
				select(NCBI_name) %>%
				unique() %>%
				pull()


print(length(TFsWithDNAbindingMotifRphaseoli))
print(TFsWithDNAbindingMotifRphaseoli)



#### How many TFs have their DNA-binding motifs conserved? #####

# Conservation ==> identity >= 40 and coverage >= 80 in both directions

print("Number of TFs with all their DNA-binding motifs conserved")

conservation_1 <- conservation %>%
						filter(motifDesc == "DNA-Binding-Region") %>%
						filter(Motif_Conservation == 1) %>%
						select(EC_locusTag) %>%
						pull() %>%
						unique()

conservation_0 <- conservation %>%
						filter(motifDesc == "DNA-Binding-Region") %>%
						filter(Motif_Conservation == 0) %>%
						select(EC_locusTag) %>%
						pull() %>%
						unique()

print(length(setdiff(conservation_1, conservation_0)))
print(setdiff(conservation_1, conservation_0))



#### What's the average identity and coverage of those motifs that are conserved? #####

print("the average identity and coverage of the conserved DNA-binding region in Ecoli is:")

ConservedDNAbindingMotif_metics <- ecoli %>%
										filter(EC_locusTag %in% setdiff(conservation_1, conservation_0)) %>%
										filter(motifDesc == "DNA-Binding-Region") %>%
										select(motifCoverage, identPercent)

print(mean(ConservedDNAbindingMotif_metics$motifCoverage))
print(mean(ConservedDNAbindingMotif_metics$identPercent))


print("the average identity and coverage of the conserved DNA-binding region in Rphaseoli is:")

ConservedDNAbindingMotif_metics_Rphaseoli <- rphaseoli %>%
										filter(EC_locusTag %in% setdiff(conservation_1, conservation_0)) %>%
										filter(motifDesc == "DNA-Binding-Region") %>%
										select(motifCoverage, identPercent)

print(mean(ConservedDNAbindingMotif_metics_Rphaseoli$identPercent))
print(mean(ConservedDNAbindingMotif_metics_Rphaseoli$motifCoverage))