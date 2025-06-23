# Name:
# section1-Figures-supp.R
#Author:
# Hernandez Montserrat
#Version
# v1
#Description
# The script allows all BLAST quality plots
# Input parameters:
#   --inpath[directory]
#   --outpath[directory]

#  Output
#   1) A set of png files with BLAST quality properties (dot plot identity vs coverage and density plot of e-value and qlen)
# 
# Example:
# Rscript --vanilla section1-Figures-supp.R
# /home/emhernan/1_BBH_TFs/BLASTresults/
# /home/emhernan/1_BBH_TFs/png/
# Rscript --vanilla section1-Figures-supp.R /home/emhernan/1_BBH_TFs/BLASTresults/ /home/emhernan/1_BBH_TFs/png/


# inpath <- "/home/emhernan/1_BBH_TFs/BLASTresults/"
# outpath <- "/home/emhernan/1_BBH_TFs/png/"


############################## Functions ###############################
###########################################################################
########################## density_plot ###################################
###########################################################################

density_plot <- function(df, x_axe, col, filler, xilim, xslim, yilim, yslim, ylab, xlab, flagcolab = FALSE, colab = "", titlelab, xsize, ysize, textsize, titlesize = 12, angleV = 45){
  
  ggobject <- df %>%
    ggplot(aes(x = x_axe)) +
    geom_density(colour = col, fill = filler)
  ggobject <- plot_settings(ggobject, xilim, xslim, yilim, yslim, ylab, xlab, flagcolab, colab, titlelab, xsize, ysize, textsize, titlesize, angleV)
  
  return(ggobject)
}

###########################################################################
######################### plot_settings ###################################
###########################################################################

plot_settings <- function(ggobj, xilim, xslim, yilim, yslim, ylab, xlab, flagcolab, colab, titlelab, xsize, ysize, textsize, titlesize, angleV){
  
  ggobj <- ggobj +
    coord_cartesian(ylim = c(yilim, yslim), xlim = c(xilim, xslim)) +
    theme_classic() +
    labs(y = ylab, 
         x = xlab,
         title = titlelab) +
    theme(plot.title = element_text(hjust = 0.5, size = titlesize, face = "italic"), 
          axis.title.x = element_text(size = xsize),
          axis.title.y = element_text(size = ysize),
          axis.text = element_text(size = textsize),
          axis.text.x = element_text(angle = angleV), 
          panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent', color=NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_colour_viridis(option="mako", direction = -1) 
  
  if(flagcolab){
    ggobj <- ggobj +
      labs(color = colab)
  }
  
  return(ggobj)
}


############################## Main program ###############################

############################ Loading libraries ############################
print("........................Loading libraries........................") 

library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(viridis)

###########################################################################
############################## Inputs #####################################
###########################################################################

arg = commandArgs(trailingOnly = T)

if (length(arg)==0) {
  stop("Must supply input and output directories", call.=FALSE)
} else if (length(arg) == 1){
  stop("Enter output directory")
}


inpath <- arg[1]
outpath <- arg[2]


print("........................Here is your input data........................")

print(inpath)
print(outpath)


print("........................Reading data........................")

ecoliTFs     <- read.table(file = paste(inpath, "TFs_ECaaq_RZaadb_blastP_b1_m8.tab", sep = ""), header =  TRUE, sep = "\t")
rphaseoliTFs <- read.table(file = paste(inpath, "TFs_RZaaq_ECaadb_blastP_b1_m8.tab", sep = ""), header =  TRUE, sep = "\t")

print("........................Processing data ........................")

ecoliEvalue <- density_plot(df = ecoliTFs, x_axe = ecoliTFs$Evalue, col = "#e2127f", filler = "#fbd1e2", 
             xilim = 0, xslim = max(ecoliTFs$Evalue), yilim = 0, yslim = 1.8e+30, ylab = "Density", xlab = "E-value", 
           titlelab = "Escherichia coli", xsize = 4.5, ysize = 4.5, textsize = 3, titlesize = 5, angleV = 45)

rphaseoliEvalue <- density_plot(df = rphaseoliTFs, x_axe = rphaseoliTFs$Evalue, col = "#e2127f", filler = "#fbd1e2", 
                                 xilim = 0, xslim =  max(rphaseoliTFs$Evalue), yilim = 0, yslim = 1.8e+30, ylab = "Density", xlab = "E-value", 
                                 titlelab = "Rhizobium phaseoli", xsize = 4.5, ysize = 4.5, textsize = 3, titlesize = 5, angleV = 45)

ecoliQlen <-  density_plot(df = ecoliTFs, x_axe = ecoliTFs$qlen, col = "#10c3a6", filler = "#93f5e5", 
                            xilim = 0, xslim = 1300, yilim = 0, yslim = 0.004, ylab = "Density", xlab = "Protein length (aa)", 
                            titlelab = "Escherichia coli", xsize = 4.5, ysize = 4.5, textsize = 3, titlesize = 5, angleV = 0)

rphaseoliQlen <-  density_plot(df = rphaseoliTFs, x_axe = rphaseoliTFs$qlen, col = "#10c3a6", filler = "#93f5e5", 
                            xilim = 0, xslim =  1300, yilim = 0, yslim = 0.004, ylab = "Density", xlab = "Protein length (aa)", 
                            titlelab = "Rhizobium phaseoli", xsize = 4.5, ysize = 4.5, textsize = 3, titlesize = 5, angleV = 0)


SuppF1 <- ggarrange(ecoliEvalue, rphaseoliEvalue, ecoliQlen, rphaseoliQlen,
                     labels = c("A", "B", "C", "D"), 
                     ncol = 2, nrow = 2)

print("........................Saving data ........................")

pngpath <- paste(outpath, "SuppF2.png",sep="")
png(pngpath, width=3*300, height=4.25*300, res=300)
SuppF1
dev.off()


