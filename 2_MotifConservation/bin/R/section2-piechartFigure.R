# Name:
#   section2-piechartFigure.R
# Author:
#   Hernandez Benitez Ericka Montserrat
# Version:
#   v1
# Description:
#   The script computes the pie chart of the motifs' conservation
# Input parameters:
#   --inpath[directory]
#   --outpath[directory]
# Output
#     1) A png file with 5 grids
#
# Rscript --vanilla section2-piechartFigure.R
# /home/emhernan/2_MotifConservation/png/
# /home/emhernan/2_MotifConservation/motifsInfo/motifsConservationIdent30.tsv
# Rscript --vanilla section2-piechartFigure.R /home/emhernan/2_MotifConservation/png/ /home/emhernan/2_MotifConservation/motifsInfo/motifsConservationIdent30.tsv

# outpath <- "/home/emhernan/2_MotifConservation/png/"
# inpath <- "/home/emhernan/2_MotifConservation/motifsInfo/motifsConservationIdent40Coverage80.tsv"


print(".................................Loading libraries.................................")


library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(png)


############################## Functions ###############################

###########################################################################
############################ pie_plot #####################################
###########################################################################

pie_plot <- function(df_tmp= df, unit = 1.5, legend_text_size = 15, legend_title_size = 20, plot_title_size = 30, vector=c("#999999", "#E69F00", "#56B4E9"), palete = FALSE){
  
  ggobject <- ggplot(df_tmp, aes(x="", y=prop, fill=Group)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() +
    theme(legend.key.size = unit(unit, 'cm'), legend.text = element_text(size=legend_text_size), legend.title = element_text(size= legend_title_size), plot.title = element_text(hjust = 0.5,  size = plot_title_size))
  
  if(palete) {
    ggobject <- ggobject +  
      scale_fill_brewer(palette="Dark2") 
  }
  
  else {
    ggobject <- ggobject +  
      scale_fill_manual(values=vector)
  }
  
  
  return(ggobject)
  
}




process_df <- function(df_tmp= df, filter_array = all, filter_num = 0, atLeast = FALSE, created = FALSE){
  
  if(atLeast){
    df_tmp <- df_tmp %>% 
      filter(Motif_Conservation == filter_num)
  }
  
  
  if(!created){
    df_tmp <- df_tmp %>% 
      filter(EC_locusTag %in% filter_array) %>%
      select(motifDesc) %>% 
      count(motifDesc) %>% 
      rename(Group = motifDesc,value = n)
    
    
  }
  df_tmp <- df_tmp %>%
    arrange(desc(Group)) %>%
    mutate(prop = value / sum(df_tmp$value) *100) %>%
    mutate(Group= paste(Group, "\n(", as.character(round(prop,2)), ")%", sep = ""))
  return(df_tmp)
}



###########################################################################
############################## Inputs #####################################
###########################################################################

arg = commandArgs(trailingOnly = T)

if (length(arg)==0) {
  stop("Must supply input and output directories", call.=FALSE)
} else if (length(arg) == 1){
  stop("Enter output directory")
}

inpath   <- arg[1]
outpath  <- arg[2]


print(".......................Reading data.......................")

TFs     <- read.table(file = inpath, header = TRUE, sep = "\t")

##### Pie chart ############

n_0   <- TFs %>% filter(Motif_Conservation == 0) %>% select(EC_locusTag) %>% unique()
n_1   <- TFs %>% filter(Motif_Conservation == 1) %>% select(EC_locusTag) %>% unique()


all     <- setdiff(n_1,n_0)$EC_locusTag
atLeast1  <- intersect(n_0,n_1)$EC_locusTag
none0s    <- setdiff(n_0,n_1)$EC_locusTag 


df <- data.frame(
  Group= c("All conserved motifs", "At least one conserved motif", "No conserved motifs"),
  value=c(length(all), length(atLeast1), length(none0s))
)

df <- process_df(df_tmp =df, atLeast = FALSE, created = TRUE)

df_All    <- process_df(df_tmp = TFs, filter_array = all, atLeast = FALSE, created = FALSE)
df_AtLeast  <- process_df(df_tmp = TFs, filter_array = atLeast1, filter_num = 1, atLeast = TRUE, created = FALSE) 
df_N0s    <- process_df(df_tmp = TFs, filter_array = atLeast1, filter_num = 0, atLeast = TRUE, created = FALSE) 
df_None0s <- process_df(df_tmp = TFs, filter_array = none0s, atLeast = FALSE, created = FALSE) 


p1 <- pie_plot(df, unit = 0.5, legend_text_size = 6.3, legend_title_size = 8, plot_title_size = 8, palete = TRUE)
p2 <- pie_plot(df_All,unit = 0.5, legend_text_size = 6.3, legend_title_size = 8, plot_title_size = 8, vector = c("#67001F", "#B2182B", "#D6604D"))
p3 <- pie_plot(df_AtLeast, unit = 0.5, legend_text_size = 6.3, legend_title_size = 8, plot_title_size = 8, vector = c( "#FDDBC7", "#67001F", "#B2182B" ,"#D6604D", "#F4A582"))
p4 <- pie_plot(df_N0s, unit = 0.5, legend_text_size = 6.3, legend_title_size = 8, plot_title_size = 8,vector = c("#B2182B", "#D6604D"))
p5 <- pie_plot(df_None0s, unit = 0.5, legend_text_size = 6.3, legend_title_size = 8, plot_title_size = 8, vector = c("#67001F","#B2182B", "#D6604D"))

Piechart <- ggpubr::ggarrange(p1,p2,p3,p4,p5,
                  labels = c("A", "B","C","D","E"),
                  ncol = 2, nrow = 3)


pngpath <- paste(outpath, "SuppF3.png",sep="")
png(pngpath, width = 6.8*300, height = 4.5*300, res = 300)
Piechart
dev.off()