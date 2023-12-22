library(tidyverse)
library(pheatmap) 

col1 <- c("nodA", "nodB", "nodC", "nodS", "nodT", "nodZ", "nolO", "nifR", "nifS")

col2 <- c(1, 1, 1, 1, 1, 1, 1, 0, 0)

col3 <- c(1, 0, 0, 0, 0, 0, 0, 1, 1)

col4 <- c(1, 1, 1, 1, 0, 1, 1, 0, 0)

col5 <- c(1, 1, 1, 1, 0, 1, 1, 0, 0)


m <- cbind(col2,col3,col4,col5)
colnames(m) <- c("bean", "maize", "milpa","maize/milpa")
rownames(m) <- col1



outpath <- "/home/emhernan/4_Evidence_Levels/png"

pngpath <- paste(outpath, "/", "SuppF4.png", sep = "")
png(pngpath, width=300*5, height=300*4, res=300, units="px", bg = "transparent")

pheatmap(m, cluster_rows = FALSE,
	     fontsize_row = 8,
         fontsize_col=8, 
         color = c("#CD5C5C", "#00008B"), 
         legend = FALSE)

dev.off()



pdfpath <- paste(outpath, "/", "SuppF4.pdf", sep = "")
pdf(pdfpath, width=5, height=4, bg = "transparent")

pheatmap(m, cluster_rows = FALSE,
	     fontsize_row = 8,
         fontsize_col=8, 
         color = c("#CD5C5C", "#00008B"), 
         legend = FALSE)

dev.off()

