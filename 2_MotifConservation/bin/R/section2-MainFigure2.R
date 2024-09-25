# Name:
#   section2-MainFigure2.R
# Author:
#   Hernandez Benitez Ericka Montserrat
# Version:
#   v1
# Description:
#   The script computes the dot plot identity vs. coverage of DNA binding motifs.
# Input parameters:
#   --inpath[directory]
#   --outpath[directory]
#   --ecoliFile[file name]
#   --rphaseoliFile[file name]
#   --color[color palette]
#   --outFile[directory]
#   --interceptT[boolean]

# Output
#     1) A png file with 2 grids

# Rscript --vanilla section2-MainFigure2.R
# /home/emhernan/2_MotifConservation/motifsInfo/
# /home/emhernan/2_MotifConservation/png/
# ECaaq_oTFs_DNAbinding_motifs_info.tsv
# RZaaq_oTFs_DNAbinding_motifs_info.tsv
# mako
# Rscript --vanilla section2-MainFigure2.R /home/emhernan/2_MotifConservation/motifsInfo/ /home/emhernan/2_MotifConservation/png/ ECaaq_oTFs_DNAbinding_motifs_info.tsv RZaaq_oTFs_DNAbinding_motifs_info.tsv mako


# inpath <- "/home/emhernan/2_MotifConservation/motifsInfo/"
# outpath <- "/home/emhernan/2_MotifConservation/png/"
# ecoliFile <- "ECaaq_oTFs_DNAbinding_motifs_info.tsv"
# rphaseoliFile <- "RZaaq_oTFs_DNAbinding_motifs_info.tsv"
# color <- "mako"


############################## Functions ###############################

###########################################################################
############################ dot_plot #####################################
###########################################################################

dot_plot <- function(df, x_axe, y_axe, xilim, xslim, yilim, yslim, ylab = "", xlab = "",textsize, scaleColor = "mako", intercept = FALSE, jitter = 1.5){
  
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
library(png)

arg = commandArgs(trailingOnly = T)

if (length(arg)==0) {
  stop("Must supply input and output directories", call.=FALSE)
} else if (length(arg) == 1){
  stop("Enter output directory")
} else if (length(arg) == 2){
  stop("Enter ecoli file name")
} else if (length(arg) == 3){
  stop("Enter rphaseoli file name")
} else if (length(arg) == 4){
  stop("Enter color palette name")
}

inpath        <- arg[1]
outpath       <- arg[2]
ecoliFile     <- arg[3]
rphaseoliFile <- arg[4]
color         <- arg[5]

print("........................Reading data........................")

ecoli     <- read.table(file = paste(inpath, ecoliFile, sep = ""), header =  FALSE, sep = "\t")
rphaseoli <- read.table(file = paste(inpath, rphaseoliFile, sep = ""), header =  FALSE, sep = "\t")


print("........................Processing data ........................")


ecoli <- ecoli %>% rename("NCBI_name" =  V1, "EC_locusTag" = V2, "motifCoverage" = V3, "identPercent" = V4, "motifDesc" = V5)
rphaseoli <- rphaseoli %>% rename("NCBI_name" =  V1, "EC_locusTag" = V2, "motifCoverage" = V3, "identPercent" = V4, "motifDesc" = V5)
  

dotplotEcoli    <- dot_plot(df = ecoli, x_axe = ecoli$motifCoverage, y_axe = ecoli$identPercent, 
                            xilim = 80, xslim = 100, yilim = 40, yslim = 80, 
                            textsize = 6, scaleColor = color, jitter = 0.5)

dotplotRphaseoli <- dot_plot(df = rphaseoli, x_axe = rphaseoli$motifCoverage, y_axe = rphaseoli$identPercent, 
                             xilim = 80, xslim = 100, yilim = 40, yslim = 80, 
                             textsize = 6, scaleColor = color, intercept = interceptT, jitter = 0.5)


print("........................Saving data ........................")

pngpath <- paste(outpath, "dotplotEcoli_motifs.png",sep="")
png(pngpath, width=1.5*300, height=2.125*300, res=300)
dotplotEcoli
dev.off()

pngpath <- paste(outpath, "dotplotRphaseoli_motifs.png",sep="")
png(pngpath, width=1.5*300, height=2.125*300, res=300)
dotplotRphaseoli
dev.off()

