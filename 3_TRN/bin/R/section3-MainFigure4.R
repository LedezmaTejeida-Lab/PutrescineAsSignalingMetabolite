# Name:
#   section3-MainFigure4.R
# Author:
#   Hernandez Benitez Ericka Montserrat
# Version:
#   v2
# Description:
#   The script computes the bar plot of TGs per each TF in the TRN. Each bar shows the number of orthologous and non-orthologous TGs
# Input parameters:
#   --networkPath[file path]
#   --outpath[directory]


# Output
#     1) A png file with a bar plot

# Rscript --vanilla section3-MainFigure4.R 
# /space24/PGC/emhernan/3_TRN/tables/MainTab1.tmp
# /space24/PGC/emhernan/3_TRN/png/
# Rscript --vanilla section3-MainFigure4.R  /space24/PGC/emhernan/3_TRN/tables/MainTab1.tmp /space24/PGC/emhernan/3_TRN/png/



# networkPath     <- "/space24/PGC/emhernan/3_TRN/tables/MainTab1.tmp"
# outpath     <- "/space24/PGC/emhernan/3_TRN/png/"



############################ Loading libraries ############################
print("........................Loading libraries........................") 

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)


arg = commandArgs(trailingOnly = T)

if (length(arg)==0) {
  stop("Must supply networkPath and orthologousPath file paths", call.=FALSE)
} else if (length(arg) == 1){
  stop("Enter outpath directory")
} 

networkPath     <- arg[1]
outpath         <- arg[2]


print("........................Here is your input data........................")

print(networkPath)
print(outpath)


print("........................Reading data........................")

network  <- read.table(file = networkPath, header =  TRUE, sep = "\t")

print("........................Processing data ........................")

df1 <- network %>%
  rename(orthologue = number_orthologous_tgs) %>%
  mutate(non_orthologue = number_of_tgs-orthologue) %>%
  select(TFname, orthologue, non_orthologue) %>%
  pivot_longer(cols = orthologue:non_orthologue,
              names_to = "type", 
              values_to = "n")
df1$type[df1$type == "non_orthologue"] <- "non-orthologue"

ggobject <- ggplot(data=df1, aes(x=TFname, y=n, fill= type)) +
  geom_bar(position = "stack", stat="identity") +
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  theme_classic() +
  labs(y = "Target Genes", 
       x = "", 
       fill = "") +
  theme(axis.title.x = element_text(size = 9.5),
        axis.title.y = element_text(size = 9.5),
        axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size = 8),
        legend.text=element_text(size=8), 
        legend.position=c(0.215, 0.95),
        legend.key.size = unit(0.45, "cm"))
print("........................Saving data ........................")

pngpath <- paste(outpath, "MainF4.png",sep="")
png(pngpath, width=300*3.25, height=300*4, res=300)
ggobject
dev.off()
