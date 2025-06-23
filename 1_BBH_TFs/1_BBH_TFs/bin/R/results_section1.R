## Main results from section 1 of article

library(dplyr)

ecoli <- read.table(file = "/home/emhernan/1_BBH_TFs/BLASTresults/TFs_ECaaq_RZaadb_blastP_b1_m8.tab", header = TRUE, sep = "\t")

rphaseoli <- read.table(file = "/home/emhernan/1_BBH_TFs/BLASTresults/TFs_RZaaq_ECaadb_blastP_b1_m8.tab", header = TRUE, sep = "\t")



print("Query length of Ecoli and Rphaseoli")

print(mean(ecoli$qlen))
print(mean(rphaseoli$qlen))

ecoli <- ecoli %>%
	mutate(qName = sub("\\|[^|]+$", "", qName)) %>%
	select(qName, sName, peri, coveragePercent) %>%
	rename("periEcoli" = peri, "coveragePercentEcoli" = coveragePercent)

rphaseoli <- rphaseoli %>%
				select(qName, sName, peri, coveragePercent) %>%
				rename("periRphaseoli" = peri, "coveragePercentRphaseoli" = coveragePercent) 


print("TFs with identity  > 50 and coverage > 50 in both directions")

TFs5050 <- ecoli %>%
	full_join(rphaseoli, by = c( "qName" =  "sName")) %>%
	filter(periEcoli > 50 & periRphaseoli > 50) %>%
	filter(coveragePercentEcoli > 50 & coveragePercentRphaseoli >50)

print(TFs5050)


print("TFs with max identity in both directions")

maxIdentity <- ecoli %>%
	full_join(rphaseoli, by = c( "qName" =  "sName")) %>%
	filter(periEcoli == max(periEcoli) | periRphaseoli == max(periRphaseoli))



print(maxIdentity)
print("TFs with max coverage in both directions")

maxCoverage <- ecoli %>%
	full_join(rphaseoli, by = c( "qName" =  "sName")) %>%
	filter(coveragePercentEcoli == max(coveragePercentEcoli) | coveragePercentRphaseoli == max(coveragePercentRphaseoli))

print(maxCoverage)