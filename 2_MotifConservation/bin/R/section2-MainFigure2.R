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
# MainF2.png
# TRUE
# Rscript --vanilla section2-MainFigure2.R /home/emhernan/2_MotifConservation/motifsInfo/ /home/emhernan/2_MotifConservation/png/ ECaaq_oTFs_DNAbinding_motifs_info.tsv RZaaq_oTFs_DNAbinding_motifs_info.tsv mako MainF2.png TRUE

# inpath <- "/home/emhernan/2_MotifConservation/motifsInfo/"
# outpath <- "/home/emhernan/2_MotifConservation/png/"
# ecoliFile <- "ECaaq_oTFs_DNAbinding_motifs_info.tsv"
# rphaseoliFile <- "RZaaq_oTFs_DNAbinding_motifs_info.tsv"
# color <- "mako"
# outFile <- "MainF2.png"
# interceptT <- FALSE


############################## Functions ###############################

###########################################################################
############################ dot_plot #####################################
###########################################################################

dot_plot <- function(df, x_axe, y_axe, xilim, xslim, yilim, yslim, ylab, xlab, flagcolab = TRUE, colab = "", titlelab, xsize, ysize, textsize, titlesize, scaleColor = "mako", intercept = FALSE, legendsize){
  
  ggobject <- df %>% 
    ggplot(aes(x = x_axe, y = y_axe, color = x_axe)) +
    geom_point() +
    geom_jitter(width = 0.4, height = 0.4) 
    ggobject <- plot_settings(ggobject, xilim, xslim, yilim, yslim, ylab, xlab, flagcolab, colab, titlelab, xsize, ysize, textsize, titlesize, scaleColor, legendsize)
  
  if(intercept){
  ggobject <- ggobject + 
              geom_vline(xintercept = 30, colour = "#000000", linetype="dashed")
  }
  return(ggobject)
  
}

###########################################################################
######################### plot_settings ###################################
###########################################################################

plot_settings <- function(ggobj, xilim, xslim, yilim, yslim, ylab, xlab, flagcolab, colab, titlelab, xsize, ysize, textsize, titlesize, scaleColor, legendsize){
  
  ggobj <- ggobj +
    coord_cartesian(ylim = c(yilim, yslim), xlim = c(xilim, xslim)) +
    theme_classic() +
    labs(y = ylab, 
         x = xlab,
         title = titlelab) +
    theme(
          legend.title=element_text(size=legendsize),
          legend.text = element_text(size = legendsize),
          plot.title = element_text(hjust = 0.5, size = titlesize, face = "italic"), 
          axis.title.x = element_text(size = xsize, hjust = 0.5, vjust = 0.5),
          axis.title.y = element_text(size = ysize, angle = 90, hjust = 0.5, vjust = 0.5),
          axis.text = element_text(size = textsize)) +
    scale_colour_viridis(option= scaleColor, direction = -1) 
  
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
  stop("Enter colo palette name")
} else if (length(arg) == 5){
  stop("Enter output file name")
} else if (length(arg) == 6){
  stop("Enter TRUE/FALSE for intercept line")
} 

inpath        <- arg[1]
outpath       <- arg[2]
ecoliFile     <- arg[3]
rphaseoliFile <- arg[4]
color         <- arg[5]
outFile       <- arg[6]
interceptT    <- arg[7]


print("........................Reading data........................")

ecoli     <- read.table(file = paste(inpath, ecoliFile, sep = ""), header =  FALSE, sep = "\t")
rphaseoli <- read.table(file = paste(inpath, rphaseoliFile, sep = ""), header =  FALSE, sep = "\t")


print("........................Processing data ........................")


ecoli <- ecoli %>% rename("NCBI_name" =  V1, "EC_locusTag" = V2, "motifCoverage" = V3, "identPercent" = V4, "motifDesc" = V5)
rphaseoli <- rphaseoli %>% rename("NCBI_name" =  V1, "EC_locusTag" = V2, "motifCoverage" = V3, "identPercent" = V4, "motifDesc" = V5)
  

dotplotEcoli    <- dot_plot(df = ecoli, x_axe = ecoli$motifCoverage, y_axe = ecoli$identPercent, 
                         xilim = 80, xslim = 100, yilim = 30, yslim = 100, ylab = "Identity (%)", 
                         xlab = "Query coverage (%)", colab = 'Query coverage (%)', 
                         titlelab = "Escherichia coli\n", xsize = 8.5, ysize = 8.5, textsize = 8.5, titlesize = 8, scaleColor = color, intercept = interceptT, legendsize = 8)

dotplotRphaseoli <- dot_plot(df = rphaseoli, x_axe = rphaseoli$motifCoverage, y_axe = rphaseoli$identPercent, 
                         xilim = 80, xslim = 100, yilim = 30, yslim = 100, ylab = "Identity (%)", 
                         xlab = "Query coverage (%)", colab = 'Query coverage (%)', 
                         titlelab = "Rhizobium phaseoli\n", xsize = 8.5, ysize = 8.5, textsize = 8.5, titlesize = 8, scaleColor = color, intercept = interceptT, legendsize = 8)


MainF2 <- ggarrange(dotplotEcoli, dotplotRphaseoli, 
                    labels = c("A", "B"), 
                    ncol = 2, nrow = 1, 
                    common.legend = TRUE, legend="bottom")


print("........................Saving data ........................")

pngpath <- paste(outpath, outFile,sep="")
png(pngpath, width=3*300, height=4.25*300, res=300)
MainF2
dev.off()
