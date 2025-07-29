## Loading libraries
library(dplyr)
library(ggvenn)
library(ggplot2)
library(tidyr)

## Loading data

bean_up <- read.table(file = "/home/emhernan/4_Evidence_Levels/output/bean_upregulated_lfc1_v2.tsv", sep = "\t", header = TRUE)
bean_down <- read.table(file = "/home/emhernan/4_Evidence_Levels/output/bean_downregulated_lfc1_v2.tsv", sep = "\t", header = TRUE)
TRN <- read.table(file = "/home/emhernan/4_Evidence_Levels/input/oTF-TGs-TRN-bgM1.tsv", sep = "\t", header = TRUE)
annnotation <- read.table(file = "/home/emhernan/4_Evidence_Levels/tables/SuppTable2.tsv", sep = "\t", header = TRUE)
all_DEA_genes <- read.table(file = "/home/emhernan/4_Evidence_Levels/output/DEA_gene_results_bean.tsv", sep = "\t", header = TRUE)
trn_tgs <- unique(unlist(strsplit(TRN$TGs, split = ",")))


## 1) Let's plot the Venn Diagramm with all upregulated genes
#============================================================================================#
list_orf_up <- list(
  expression = unique(c(bean_up$Locus)),
  trn = trn_tgs,
  function_sim =  unique(annnotation$locusTag)
  
)

png(filename = "/home/emhernan/4_Evidence_Levels/png/upregulated_evidencelevel_diagramm.png", width=7*300, height=7*300, res=300)
names(list_orf_up) <- c("Expression", "Transcriptional Regualtion", "Function")
ggvenn(list_orf_up, fill_color = c("#E44F2C","#21908CFF", "#7570B3"), stroke_size = 0.5, set_name_size = 7, text_size = 8, show_percentage = FALSE)

dev.off()

## 2) Let's plot the Venn Diagramm with all downregulated genes
#============================================================================================#

list_orf_down <- list(
  expression = unique(c(bean_down$Locus)),
  trn = trn_tgs,
  function_sim =  unique(annnotation$locusTag)

  )

png(filename = "/home/emhernan/4_Evidence_Levels/png/downregulated_evidencelevel_diagramm.png", width=7*300, height=7*300, res=300)

names(list_orf_down) <- c("Expression", "Transcriptional Regualtion", "Function")
ggvenn(list_orf_down, fill_color = c("#E44F2C","#21908CFF", "#7570B3"), stroke_size = 0.5, set_name_size = 7, text_size = 8, show_percentage = FALSE)

dev.off()

## 3) Let's plot the Venn Diagramm with all DEA genes
#============================================================================================#
list_orf_allDE <- list(
  expression =  unique(c(bean_down$Locus, bean_up$Locus))[!unique(c(bean_down$Locus, bean_up$Locus)) == ""],
  trn = trn_tgs,
  function_sim =  unique(annnotation$locusTag)
  
)

png(filename = "/home/emhernan/4_Evidence_Levels/png/all_de_regulated_diagramm.png", width=7*300, height=7*300, res=300)
names(list_orf_allDE) <- c("Expression", "Transcriptional Regulation", "Function")
ggvenn(list_orf_allDE, fill_color = c("#E44F2C","#21908CFF", "#7570B3"), stroke_size = 0.5, set_name_size = 7, text_size = 8, show_percentage = FALSE)
dev.off()


#========================================================================================#
########################### Metrics and tables for the article  ##########################
#========================================================================================#


#========================================================================================#
######################### Literature and differential expression #########################
#========================================================================================#

print("Literature and differential expression intersection - it contemplates genes in the trn intersection")
print("These genes with known symbiosis function are down regulated")
df3 <- annnotation %>%
  filter(locusTag %in% intersect(list_orf_down$Expression, list_orf_down$Function)) %>%
  arrange(geneName)
print(df3)


print("These genes with known symbiosis function are up regulated ")
df4 <- annnotation %>%
  filter(locusTag %in% intersect(list_orf_up$Expression, list_orf_up$Function)) %>%
  arrange(geneName)  
print(df4)

df5 <- rbind(df3, df4) %>% 
            left_join(all_DEA_genes, by = c ("locusTag" = "Locus")) %>%
            select(geneName, locusTag, logFC, adj.PVal, FDR) %>%
            arrange(geneName)

write.table(x = df5, file = "/home/emhernan/4_Evidence_Levels/output/nodGenesDEA_info.tsv", sep =  "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



#========================================================================================#
######################### Literature,  differential expression, and trn #########################
#========================================================================================#


names(list_orf_allDE) <- c("expression", "trn", "function_sim")
inter <- intersect(intersect(list_orf_allDE$expression, list_orf_allDE$trn), list_orf_allDE$function_sim)

TRN_formatted <- TRN %>%
    select(TFname, TGs) %>%
    separate_rows(TGs, sep = ",") %>%
    distinct()

tf_tg_lit_dea_trn_df <- annnotation %>%
  filter(locusTag %in% inter) %>%
  arrange(geneName) %>%
  left_join(TRN_formatted, by = c("locusTag" = "TGs")) %>%
  select(TFname, geneName) %>%
  group_by(TFname) %>%
  summarise(TGs = paste(geneName, collapse = ","), .groups = "drop")

n <- annnotation %>% filter(locusTag %in% inter) %>% select(geneName) %>% pull() %>% length()
print(paste("There are:", as.character(n), "genes in the intersection of function, expression and trnr"))
print("The genes are:")
print(tf_tg_lit_dea_trn_df)
#========================================================================================#
################################### Literature and DEA ###################################
#========================================================================================#
lit_dea <- setdiff(intersect(list_orf_allDE$expression, list_orf_allDE$function_sim), inter)
lit_dea_df <- annnotation %>%
              filter(locusTag %in% lit_dea) %>%
              arrange(geneName) 
print(paste("There are", as.character(nrow(lit_dea_df)), "genes in the function and Esxpression intersection"), sep = " ")
print("The genes are:")
print(lit_dea_df)

#========================================================================================#
################################### Literature and trn ###################################
#========================================================================================#

lit_tnr <- setdiff(intersect(list_orf_allDE$trn, list_orf_allDE$function_sim), inter)
lit_trn_df <- annnotation %>%
              filter(locusTag %in% lit_tnr) %>%
              arrange(geneName) 
print(paste("There are", as.character(nrow(lit_trn_df)), "genes in the function and trn intersection"), sep = " ")
print("The genes are:")
print(lit_trn_df)
print("And are regulated by:")

tf_tg_lit_trn_df <- lit_trn_df %>%
      left_join(TRN_formatted, by = c("locusTag" = "TGs")) %>%
      select(TFname, geneName) %>%
      group_by(TFname) %>%
      summarise(TGs = paste(geneName, collapse = ","), .groups = "drop")
print(tf_tg_lit_trn_df)

write.table(tf_tg_lit_trn_df, file = "/home/emhernan/4_Evidence_Levels/tables/SuppTable3.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
