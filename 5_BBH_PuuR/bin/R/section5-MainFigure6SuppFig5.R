# Name:
#   section5-MainFigure6SuppFig5.R
# Author:
#   Hernandez Benitez Ericka Montserrat
# Version:
#   v1
# Description:
#   The script computes the dot plot identity vs. coverage of PuuR BBH and the E-value density.
# Input parameters:
#   --inpath[directory]
#   --outpath[directory]

# Output
#     1) Two png files 

# Rscript --vanilla section5-MainFigure6SuppFig5.R
# /home/emhernan/5_BBH_PuuR/output/
# /home/emhernan/5_BBH_PuuR/png/
# Rscript --vanilla section5-MainFigure6SuppFig5.R /home/emhernan/5_BBH_PuuR/output/ /home/emhernan/5_BBH_PuuR/png/

# inpath <- "/home/emhernan/5_BBH_PuuR/output/"
# outpath <- "/home/emhernan/5_BBH_PuuR/png/"

############################## Functions ###############################

###########################################################################
############################ dot_plot #####################################
###########################################################################

dot_plot <- function(df, x_axe, y_axe, xilim, xslim, yilim, yslim, ylab, xlab, flagcolab = TRUE, colab = "", titlelab, xsize, ysize, textsize, titlesize, legendsize = 7, angleV = 0){
  
  ggobject <- df %>% 
    ggplot(aes(x = x_axe, y = y_axe, color = x_axe)) +
    geom_point(position=position_jitter(h=0.7, w=0.7), size = 5) 
  ggobject <- plot_settings(ggobject, xilim, xslim, yilim, yslim, ylab, xlab, flagcolab, colab, titlelab, xsize, ysize, textsize, titlesize, legendsize, angleV)
  
  return(ggobject)
  
}


###########################################################################
########################## density_plot ###################################
###########################################################################

density_plot <- function(df, x_axe, col, filler, xilim, xslim, yilim, yslim, ylab, xlab, flagcolab = FALSE, colab = "", titlelab, xsize, ysize, textsize, titlesize, legendsize = 7, angleV = 0){
  
  ggobject <- df %>%
    ggplot(aes(x = x_axe)) +
    geom_density(colour = col, fill = filler)
  ggobject <- plot_settings(ggobject, xilim, xslim, yilim, yslim, ylab, xlab, flagcolab, colab, titlelab, xsize, ysize, textsize, titlesize, legendsize, angleV)
  
  return(ggobject)
}


###########################################################################
######################### plot_settings ###################################
###########################################################################

plot_settings <- function(ggobj, xilim, xslim, yilim, yslim, ylab, xlab, flagcolab, colab, titlelab, xsize, ysize, textsize, titlesize, legendsize, angleV){
  
  ggobj <- ggobj +
    coord_cartesian(ylim = c(yilim, yslim), xlim = c(xilim, xslim)) +
    theme_classic() +
    labs(y = ylab, 
         x = xlab,
         title = titlelab) +
    theme(legend.title=element_text(size=legendsize),
          legend.text = element_text(size = legendsize),
          plot.title = element_text(hjust = 0.5, size = titlesize, face = "italic"), 
          axis.title.x = element_text(size = xsize),
          axis.title.y = element_text(size = ysize),
          axis.text.y = element_text(size = textsize),
          axis.text.x = element_text(angle = angleV, hjust = 1, size = textsize)) +
    scale_colour_viridis(option="inferno", direction = -1) 
  
  if(flagcolab){
    ggobj <- ggobj +
      labs(color = colab)
  }
  
  return(ggobj)
}


print("........................Here is your input data........................")



############################## Main program ###############################

############################ Loading libraries ############################

print("........................Loading libraries........................") 

library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(viridis)



arg = commandArgs(trailingOnly = T)

if (length(arg)==0) {
  stop("Must supply input and output directories", call.=FALSE)
} else if (length(arg) == 1){
  stop("Enter output directory")
} 

inpath        <- arg[1]
outpath       <- arg[2]

print("........................Reading data........................")

rhizobium  <- read.table(file = paste(inpath, "Rpaaq_Rhaadb_PuuROrthologous.tsv", sep = ""), header =  FALSE, sep = "\t")
rhizobia <- read.table(file = paste(inpath, "Rhaaq_Rpaadb_PuuROrthologous.tsv", sep = ""), header =  FALSE, sep = "\t")

colnames(rhizobium) <- c("qName", "sName", "peri", "alilen", "numMM", "nnGP", "qSS", "qSE", "sSS", "sSE", "Evalue", "bitScore")
colnames(rhizobia) <- c("qName", "sName", "peri", "alilen", "numMM", "nnGP", "qSS", "qSE", "sSS", "sSE", "Evalue", "bitScore")

rhizobium$qcovs = 100
rhizobia$qcovs = 100

print("........................Processing data ........................")
dotplotRhizobium <- dot_plot(df = rhizobium, x_axe = rhizobium$peri, y_axe = rhizobium$qcovs, 
         xilim = 75, xslim = 100, yilim = 75, yslim = 100, ylab = "Query coverage (%)", 
         xlab = "Identity (%)", colab = "Identity (%)", 
         titlelab = "Rhizobium phaseoli\n", xsize = 8.5, ysize = 8.5, textsize = 8.5, titlesize = 8.5, legendsize = 8, angleV = 45)  

dotplotRhizobia <- dot_plot(df = rhizobia, x_axe = rhizobia$peri, y_axe = rhizobia$qcovs, 
                             xilim = 75, xslim = 100, yilim = 75, yslim = 100, ylab = "Query coverage (%)", 
                             xlab = "Identity (%)", colab = "Identity (%)", 
                             titlelab = "Hyphomicrobiales\n", xsize = 8.5, ysize = 8.5, textsize = 8.5, titlesize = 8.5, legendsize = 8, angleV = 45)  


rhizobiumEvalue <- density_plot(df = rhizobium, x_axe = rhizobium$Evalue, col = "#e2127f", filler = "#fbd1e2", 
                            xilim = 0, xslim =  max(rhizobium$Evalue), yilim = 0, yslim = 1.8e+104, ylab = "Density", xlab = "E-value", 
                            titlelab = "Rhizobium phaseoli\n", xsize =7.5, ysize = 7.5, textsize = 7.5, titlesize = 6, legendsize = 8, angleV = 45)

rhizobiaEvalue <- density_plot(df = rhizobia, x_axe = rhizobia$Evalue, col = "#e2127f", filler = "#fbd1e2", 
                                xilim = 0, xslim = max(rhizobia$Evalue), yilim = 0, yslim = 1.8e+104, ylab = "Density", xlab = "E-value", 
                                titlelab = "Hyphomicrobiales\n", xsize = 7.5, ysize = 7.5, textsize = 7.5, titlesize = 6, legendsize = 8, angleV = 45)


print("........................Saving data ........................")
MainF6 <- ggarrange(dotplotRhizobium, dotplotRhizobia, 
                    labels = c("A", "B"), 
                    ncol = 2, nrow = 1, 
                    common.legend = TRUE, legend="bottom")

SuppF5 <- ggarrange(rhizobiumEvalue, rhizobiaEvalue,
                    labels = c("A", "B"), 
                    ncol = 1, nrow = 2,
                    common.legend = FALSE)


pngpath <- paste(outpath, "MainF6.png",sep="")
png(pngpath, width=3*300, height=4.25*300, res=300)
MainF6
dev.off()

pngpath <- paste(outpath, "SuppF5.png",sep="")
png(pngpath, width=3*300, height=4.25*300, res=300)
SuppF5
dev.off()
