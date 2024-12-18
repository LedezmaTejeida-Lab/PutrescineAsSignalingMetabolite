## Loading libraries
library(dplyr)
library(ggvenn)
library(ggplot2)


## Loading data

bean_up <- read.table(file = "/home/emhernan/4_Evidence_Levels/output/bean_upregulated_lfc1_v2.tsv", sep = "\t", header = TRUE)
bean_down <- read.table(file = "/home/emhernan/4_Evidence_Levels/output/bean_downregulated_lfc1_v2.tsv", sep = "\t", header = TRUE)

TRN <- read.table(file = "/home/emhernan/4_Evidence_Levels/input/oTF-TG-TRN-10-5-bgM1.tsv", sep = "\t", header = TRUE)
annnotation <- read.table(file = "/home/emhernan/4_Evidence_Levels/input/SuppTable2.tsv", sep = "\t", header = TRUE)

all_DEA_genes <- read.table(file = "/home/emhernan/4_Evidence_Levels/output/DEA_gene_results_bean.tsv", sep = "\t", header = TRUE)


TRN <- TRN %>% filter(TF_name != "FNR")

## Let's plot the Venn Diagramm with upregulated genes (all genes)
#============================================================================================#
list_orf_up <- list(
  expression = unique(c(bean_up$Locus)),
  trn = unique(TRN$TG_oldRZlocusTag),
  function_sim =  unique(annnotation$locusTag)
  
)

png(filename = "/home/emhernan/4_Evidence_Levels/png/upregulated_evidencelevel_diagramm.png", width=7*300, height=7*300, res=300)
names(list_orf_up) <- c("Expression", "Transcriptional Regualtion", "Function")
ggvenn(list_orf_up, fill_color = c("#E44F2C","#21908CFF", "#7570B3"), stroke_size = 0.5, set_name_size = 7, text_size = 8, show_percentage = FALSE)

dev.off()


print("The genes we are missing from expression and function intersection
      are nodT (RPHASCH2410_CH16550) and nifSch (RPHASCH2410_CH21455)")  #######Checar esto!!!!!


#============================================================================================#

list_orf_down <- list(
  expression = unique(c(bean_down$Locus)),
  trn = unique(TRN$TG_oldRZlocusTag),
  function_sim =  unique(annnotation$locusTag)

  )

png(filename = "/home/emhernan/4_Evidence_Levels/png/downregulated_evidencelevel_diagramm.png", width=7*300, height=7*300, res=300)

names(list_orf_down) <- c("Expression", "Transcriptional Regualtion", "Function")
ggvenn(list_orf_down, fill_color = c("#E44F2C","#21908CFF", "#7570B3"), stroke_size = 0.5, set_name_size = 7, text_size = 8, show_percentage = FALSE)

dev.off()

print("The gene we identified from downregualted analysis in the intersection is fixA (RPHASCH2410_PC00145)")
#============================================================================================#
list_orf_allDE <- list(
  expression =  unique(c(bean_down$Locus, bean_up$Locus))[!unique(c(bean_down$Locus, bean_up$Locus)) == ""],
  trn = unique(TRN$TG_oldRZlocusTag),
  function_sim =  unique(annnotation$locusTag)
  
)


png(filename = "/home/emhernan/4_Evidence_Levels/png/all_de_regulated_diagramm.png", width=7*300, height=7*300, res=300)

names(list_orf_allDE) <- c("Expression", "Transcriptional Regualtion", "Function")
ggvenn(list_orf_allDE, fill_color = c("#E44F2C","#21908CFF", "#7570B3"), stroke_size = 0.5, set_name_size = 7, text_size = 8, show_percentage = FALSE)


dev.off()


#========================================================================================#

intersect(intersect(list_orf_allDE$expression, list_orf_allDE$trn), list_orf_allDE$function_sim)
# RPHASCH2410_PC00145 fixA
# RPHASCH2410_PC00450 nodI-


######################### Literature and differential expression #########################
print("Literature and differential expression intersection")
print("These genes are down regulated")
df3 <- annnotation %>%
  filter(locusTag %in% intersect(list_orf_down$Expression, list_orf_down$Function))
print(df3)


print("These genes are up regulated ")
df4 <- annnotation %>%
  filter(locusTag %in% intersect(list_orf_up$Expression, list_orf_up$Function))
  
print(df4)

df5 <- df4 %>% 
            left_join(all_DEA_genes, by = c ("locusTag" = "Locus")) %>%
            select(geneName, locusTag, logFC, adj.PVal, FDR)
write.table(x = df5, file = "/home/emhernan/4_Evidence_Levels/output/nodGenesDEA_info.tsv", sep =  "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


######################### Literature and trn #########################
print("Literature and trn intersection")

annnotation %>%
  filter(locusTag %in% intersect(list_orf_allDE$`Transcriptional Regualtion`, list_orf_allDE$Function)) %>% select(geneName) %>% pull() %>% sort(decreasing = TRUE)

df6 <- annnotation %>%
  filter(locusTag %in% intersect(list_orf_allDE$`Transcriptional Regualtion`, list_orf_allDE$Function)) %>%
  filter(geneName != "nodI" & geneName != "fixA") %>%
  left_join(TRN %>% select(TF_name,TG_oldRZlocusTag) %>% unique(), by = c("locusTag"="TG_oldRZlocusTag")) %>%
  select(TF_name, geneName) %>%
  group_by(TF_name) %>%
  summarise(geneNames = paste(geneName, collapse = ", "), .groups = 'drop') %>%
  group_by(geneNames) %>%
  summarise(TF_name = paste(TF_name, collapse = ", "), .groups = 'drop') %>%
  select(TF_name, geneNames) %>%
  arrange(desc(geneNames))

write.table(x = df6, file = "/home/emhernan/4_Evidence_Levels/output/SuppTable4.tsv", sep =  "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)






# Remember to exclude nodI and fixA 