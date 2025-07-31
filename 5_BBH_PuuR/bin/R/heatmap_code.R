#Plot similarities between PuuR orthologues in different rhizobial species.
library(pheatmap)
data <- read.table(file="/home/emhernan/5_BBH_PuuR/docs/heatmap/oPuuR_rhizobiales.txt",
                   header=FALSE, row.names=NULL, sep="\t")
colnames(data) <- c("query", "subject", "identity", "length", "mismatches", "gap.opens", "q.start", "q.end", "s.start", "s.end", "evalue", "bit.score", "positives")
data.m <- matrix(ncol=9, nrow=9)
colnames(data.m) <- unique(data$query)
rownames(data.m) <- unique(data$query)
for(i in colnames(data.m)){
  for(j in rownames(data.m)){
    data.m[i,j] <- data[which(data$query==i & data$subject==j),"identity"]
  }
}
#Plot
labels_row = c("E. coli","B. australiense", "B. diazoefficiens", "B. japonicum", "S. fredii","S. meliloti","R. leguminosarum","R. etli", "R. phaseoli")
labels_col = c("E. coli","B. australiense", "B. diazoefficiens","B. japonicum", "S. fredii","S. meliloti","R. leguminosarum","R. etli", "R. phaseoli")
newrows <- lapply(labels_row,function(x) bquote(italic(.(x))))
newcols <- lapply(labels_col,function(x) bquote(italic(.(x))))


pngpath <- "/home/emhernan/5_BBH_PuuR/png/MainF3.png"
png(pngpath, width=300*5.5, height=300*4.5, res=300, units="px")
pheatmap(data.m, labels_row=as.expression(newrows),  
         labels_col=as.expression(newcols), display_numbers = TRUE,
         number_color = "black", 
         fontsize_number = 8, cluster_rows = F, cluster_cols = F)


dev.off()