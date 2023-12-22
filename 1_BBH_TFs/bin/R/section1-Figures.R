# Name:
# section1_Figures.R
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
# Rscript --vanilla section1_Figures.R
# /home/emhernan/1_BBH_TFs/BLASTresults/
# /home/emhernan/1_BBH_TFs/png/
# Rscript --vanilla section1_Figures.R /home/emhernan/1_BBH_TFs/BLASTresults/ /home/emhernan/1_BBH_TFs/png/


# inpath <- "/home/emhernan/1_BBH_TFs/BLASTresults/"
# outpath <- "/home/emhernan/1_BBH_TFs/png/"


############################## Functions ###############################

###########################################################################
############################ dot_plot #####################################
###########################################################################


dot_plot <- function(df, x_axe, y_axe, xilim, xslim, yilim, yslim, ylab, xlab, flagcolab = TRUE, colab = "", titlelab, xsize, ysize, textsize, titlesize = 12){
  
  ggobject <- df %>% 
    ggplot(aes(x = x_axe, y = y_axe, color = x_axe)) +
    geom_point(position=position_jitter(h=0.1, w=0.1))
  ggobject <- plot_settings(ggobject, xilim, xslim, yilim, yslim, ylab, xlab, flagcolab, colab, titlelab, xsize, ysize, textsize, titlesize)
  
  return(ggobject)
  
}

###########################################################################
########################## density_plot ###################################
###########################################################################

density_plot <- function(df, x_axe, col, filler, xilim, xslim, yilim, yslim, ylab, xlab, flagcolab = FALSE, colab = "", titlelab, xsize, ysize, textsize, titlesize = 12){
  
  ggobject <- df %>%
    ggplot(aes(x = x_axe)) +
    geom_density(colour = col, fill = filler)
  ggobject <- plot_settings(ggobject, xilim, xslim, yilim, yslim, ylab, xlab, flagcolab, colab, titlelab, xsize, ysize, textsize, titlesize)
  
  return(ggobject)
}

###########################################################################
######################### plot_settings ###################################
###########################################################################

plot_settings <- function(ggobj, xilim, xslim, yilim, yslim, ylab, xlab, flagcolab, colab, titlelab, xsize, ysize, textsize, titlesize){
  
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

ecoli     <- read.table(file = paste(inpath, "TFs_ECaaq_RZaadb_blastP_b1_m8.tab", sep = ""), header =  TRUE, sep = "\t")
rphaseoli <- read.table(file = paste(inpath, "TFs_RZaaq_ECaadb_blastP_b1_m8.tab", sep = ""), header =  TRUE, sep = "\t")

print("........................Processing data ........................")

dotplotEcoli <- dot_plot(df = ecoli, x_axe = ecoli$peri, y_axe = ecoli$coveragePercent, 
                      xilim = 0, xslim = 100, yilim = 0, yslim = 100, ylab = "Query Coverage (%)", 
                       xlab = "Identity (%)", colab = 'Identity (%)', 
                       titlelab = "Escherichia coli\n", xsize = 10, ysize = 10, textsize = 8, titlesize = 12)

dotplotRphaseoli <- dot_plot(df = rphaseoli, x_axe = rphaseoli$peri, y_axe = rphaseoli$coveragePercent, 
                          xilim = 0, xslim = 100, yilim = 0, yslim = 100, ylab = "Query Coverage (%)", 
                          xlab = "Identity (%)", colab = 'Identity (%)', 
                          titlelab = "Rhizobium phaseoli\n", xsize = 10, ysize = 10, textsize = 8, titlesize = 12)


ecoliEvalue <- density_plot(df = ecoli, x_axe = ecoli$Evalue, col = "#e2127f", filler = "#fbd1e2", 
             xilim = 0, xslim =  max(ecoli$Evalue), yilim = 0, yslim = 3.5e+23, ylab = "Density", xlab = "E-value", 
           titlelab = "Escherichia coli", xsize = 8, ysize = 8, textsize = 6, titlesize = 9)

rphaseoliEvalue <- density_plot(df = rphaseoli, x_axe = rphaseoli$Evalue, col = "#e2127f", filler = "#fbd1e2", 
                                 xilim = 0, xslim =  max(rphaseoli$Evalue), yilim = 0, yslim = 14e+23, ylab = "Density", xlab = "E-value", 
                                 titlelab = "Rhizobium phaseoli", xsize = 8, ysize = 8, textsize = 6, titlesize = 9)

ecoliQlen <-  density_plot(df = ecoli, x_axe = ecoli$qlen, col = "#10c3a6", filler = "#93f5e5", 
                            xilim = 0, xslim = 1300, yilim = 0, yslim = 0.004, ylab = "Density", xlab = "Query length (bp)", 
                            titlelab = "Escherichia coli", xsize = 8, ysize = 8, textsize = 6, titlesize = 9)

rphaseoliQlen <-  density_plot(df = rphaseoli, x_axe = rphaseoli$qlen, col = "#10c3a6", filler = "#93f5e5", 
                            xilim = 0, xslim =  1300, yilim = 0, yslim = 0.004, ylab = "Density", xlab = "Query length (bp)", 
                            titlelab = "Rhizobium phaseoli", xsize = 8, ysize = 8, textsize = 6, titlesize = 9)


MainF1 <- ggarrange(dotplotEcoli, dotplotRphaseoli, 
                          labels = c("A", "B"), 
                          ncol = 2, nrow = 1, 
                          common.legend = TRUE, legend="bottom") +
		  theme(legend.background = element_rect(fill='transparent'),
          legend.box.background = element_rect(fill='transparent'))

SuppF1 <- ggarrange(ecoliEvalue, rphaseoliEvalue, ecoliQlen, rphaseoliQlen,
                     labels = c("A", "B", "C", "D"), 
                     ncol = 2, nrow = 2)

print("........................Saving data ........................")

pngpath <- paste(outpath, "MainF1.png",sep="")
png(pngpath, width=5.59, height=6.2, res=300, units="in", bg = "transparent")
MainF1
dev.off()

pngpath <- paste(outpath, "SuppF2.png",sep="")
png(pngpath, width=4.5, height=5.5, res=300, units="in", bg = "transparent")
SuppF1
dev.off()

pdfpath <- paste(outpath, "MainF1.pdf",sep="")
pdf(pdfpath, width=5.59, height=6.2, pointsize=300, bg = "transparent")
MainF1
dev.off()

pdfpath <- paste(outpath, "SuppF2.pdf",sep="")
pdf(pdfpath, width=4.5, height=5.5, pointsize=300, bg = "transparent")
SuppF1
dev.off()

