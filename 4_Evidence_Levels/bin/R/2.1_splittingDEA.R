## This script splits the DEA table into 
## Bean, Maize and Milpa up and down regulated DE genes' tables
## The scripts assumes the DEA R project


## Loading libraries
library(dplyr)

## Loading data
DEA_genes <- read.table(file = "/home/emhernan/4_Evidence_Levels/output/DEA_gene_results_bean.tsv", header = TRUE, sep = "\t")
head(DEA_genes)


## Splitting data into upregulated and downregulated

bean_up <- DEA_genes %>%
              select(Geneid, Locus, ID, logFC) %>%
              filter(logFC > 1)

bean_down <- DEA_genes %>%
                select(Geneid, Locus, ID, logFC) %>%
                filter(logFC < -1)

print(paste("There are", as.character(length(unique(DEA_genes$Geneid))), "differentially expressed genes"))

print(paste("There are", as.character(length(unique(DEA_genes$Geneid))/6539 *100), "% differentially expressed"))

print(paste("There are", as.character(length(unique(bean_up$Geneid))), "overexpressed genes"))

print(paste("There are", as.character(length(unique(bean_down$Geneid))), "underexpressed genes"))




###########################################################

## Saving files

write.table(x = bean_up, file = "/home/emhernan/4_Evidence_Levels/output/bean_upregulated_lfc1_v2.tsv", quote =  FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = bean_down, file = "/home/emhernan/4_Evidence_Levels/output/bean_downregulated_lfc1_v2.tsv", quote =  FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

