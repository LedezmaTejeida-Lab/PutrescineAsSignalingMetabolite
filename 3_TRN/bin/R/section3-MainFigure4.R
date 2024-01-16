# Name:
#   section3-MainFigure4.R
# Author:
#   Hernandez Benitez Ericka Montserrat
# Version:
#   v1
# Description:
#   The script computes the bar plot of TGs per each TF in the TRN. Each bar shows the number of orthologous and non-orthologous TGs
# Input parameters:
#   --networkPath[file path]
#   --orthologousPath[file path]
#   --outpath[directory]


# Output
#     1) A png file with a bar plot

# Rscript --vanilla section3-MainFigure4.R 
# /space24/PGC/emhernan/3_TRN/output/oTF-TG-TRN-10-5-bgM1.tsv
# /space24/PGC/emhernan/3_TRN/OrthologousTFInfo/Orthologous_table.tsv
# /space24/PGC/emhernan/3_TRN/png/
# Rscript --vanilla section3-MainFigure4.R  /space24/PGC/emhernan/3_TRN/output/oTF-TG-TRN-10-5-bgM1.tsv /space24/PGC/emhernan/3_TRN/OrthologousTFInfo/Orthologous_table.tsv /space24/PGC/emhernan/3_TRN/png/



# networkPath     <- "/space24/PGC/emhernan/3_TRN/output/oTF-TG-TRN-10-5-bgM1.tsv"
# orthologousPath <- "/space24/PGC/emhernan/3_TRN/OrthologousTFInfo/Orthologous_table.tsv"
# outpath     <- "/space24/PGC/emhernan/3_TRN/png/"



############################ Loading libraries ############################
print("........................Loading libraries........................") 

library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)


arg = commandArgs(trailingOnly = T)

if (length(arg)==0) {
  stop("Must supply networkPath and orthologousPath file paths", call.=FALSE)
} else if (length(arg) == 1){
  stop("Enter orthologousPath file path")
} else if (length(arg) == 2){
  stop("Enter outpath directory")
} 

networkPath     <- arg[1]
orthologousPath <- arg[2]
outpath         <- arg[3]


print("........................Here is your input data........................")

print(networkPath)
print(orthologousPath)
print(outpath)


print("........................Reading data........................")

network         <- read.table(file = networkPath, header =  TRUE, sep = "\t")
orthologous     <- read.table(file = orthologousPath, header =  TRUE, sep = "\t")

print("........................Processing data ........................")


df1 <- network %>% 
  select(TF_name, TF_locusTag, TG_oldRZlocusTag) %>%
  mutate(OrthologusTG = ifelse(TG_oldRZlocusTag %in% orthologous$RZ_locusTag_old,"Ortholog","non-Ortholog")) %>%
  count(TF_name, OrthologusTG)


ggobject <- ggplot(data=df1, aes(x=TF_name, y=n, fill= OrthologusTG)) +
  geom_bar(position = "stack", stat="identity") +
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  theme_classic() +
  labs(y = "Target Genes", 
       x = "Transcription Factor", 
       fill = "") +
  theme(axis.title.x = element_text(size = 9.5),
        axis.title.y = element_text(size = 9.5),
        axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size = 8),
        legend.text=element_text(size=9.5), 
        legend.position="top")
print("........................Saving data ........................")

pngpath <- paste(outpath, "MainF4.png",sep="")
png(pngpath, width=300*3.25, height=300*4, res=300)
ggobject
dev.off()
