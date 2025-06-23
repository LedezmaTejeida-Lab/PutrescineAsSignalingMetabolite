# Name:
# section1-Figures-main.R
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
# Rscript --vanilla section1-Figures-main.R
# /home/emhernan/1_BBH_TFs/BLASTresults/
# /home/emhernan/1_BBH_TFs/png/
# Rscript --vanilla section1-Figures-main.R /home/emhernan/1_BBH_TFs/BLASTresults/ /home/emhernan/1_BBH_TFs/png/


# inpath <- "/home/emhernan/1_BBH_TFs/BLASTresults/"
# outpath <- "/home/emhernan/1_BBH_TFs/png/"


############################## Functions ###############################

###########################################################################
############################ dot_plot #####################################
###########################################################################

dot_plot <- function(df, x_axe, y_axe, xilim, xslim, yilim, yslim, ylab = "", xlab = "",textsize, scaleColor = "mako", jitter = 1.5){
  
  ggobject <- df %>% 
    ggplot(aes(x = x_axe, y = y_axe, color = x_axe)) +
    geom_point(size = 0.7, shape = 19, position = position_jitter(width=jitter, height=jitter))
     
  ggobject <- plot_settings(ggobject, xilim, xslim, yilim, yslim, ylab, xlab, textsize, scaleColor)

}  

###########################################################################
######################### plot_settings ###################################
###########################################################################

plot_settings <- function(ggobj, xilim, xslim, yilim, yslim, ylab, xlab, textsize, scaleColor){
  ggobj <- ggobj +
    coord_fixed(ratio = 1, ylim = c(yilim, yslim), xlim = c(xilim, xslim), expand = TRUE) +
    scale_x_continuous(breaks = seq(xilim, xslim, by = 10)) +
    theme_classic() +
    labs(
      y = ylab, 
      x = xlab) +
    theme(
      axis.text = element_text(size = textsize), 
      axis.text.x = element_text(angle = 0),
      legend.position="none") +
    scale_colour_viridis(option= scaleColor, direction = -1) 
  
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

dotplotEcoliTFs    <- dot_plot(df = ecoliTFs, x_axe = ecoliTFs$coveragePercent, y_axe = ecoliTFs$peri, 
                            xilim = 80, xslim = 100, yilim = 30, yslim = 80, 
                            textsize = 6, scaleColor = "mako", jitter =0.5)

dotplotRphaseoliTFs <- dot_plot(df = rphaseoliTFs, x_axe = rphaseoliTFs$coveragePercent, y_axe = rphaseoliTFs$peri, 
                             xilim = 80, xslim = 100, yilim = 30, yslim = 80, 
                             textsize = 6, scaleColor = "mako", jitter = 0.5)


print("........................Saving data ........................")

pngpath <- paste(outpath, "dotplotEcoli_TFs.png",sep="")
png(pngpath, width=1.5*300, height=2.125*300, res=300)
dotplotEcoliTFs
dev.off()

pngpath <- paste(outpath, "dotplotRphaseoli_TFs.png",sep="")
png(pngpath, width=1.5*300, height=2.125*300, res=300)
dotplotRphaseoliTFs
dev.off()


