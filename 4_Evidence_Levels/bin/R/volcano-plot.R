data <- read.delim(file="/home/emhernan/4_Evidence_Levels/input/nodGenes-function-and-DEA.tsv",header=TRUE,row.names=NULL)
data <- cbind(data,-log10(data$adjPVal))
colnames(data)[7] <- "log10Pval"
DEG <- rep(FALSE,dim(data)[1])
data <- cbind(data,DEG)
data[which(data$adjPVal <= 0.05 & data$logFC > 1),"DEG"] <- TRUE
#label <- rep(NA,dim(data)[1])
#data <- cbind(data,label)
#data[which(data$adjPVal <= 0.05 & data$logFC > 1),"label"] <- data[which(data$adjPVal <= 0.05 & data$logFC > 1),"geneName"]
library(ggplot2)
library(ggrepel)

pngpath <- "/home/emhernan/4_Evidence_Levels/png/SuppF6.png"
png(pngpath, width=300*5, height=300*4, res=300, units="px")

ggplot(data=data,aes(x=logFC,y=log10Pval,label=geneName, color=DEG)) + geom_point() + 
  scale_x_continuous(name ="log2 Fold Change", limits=c(-5,5)) +
  theme_classic() + geom_text_repel() +  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  scale_color_manual(values=c("grey", "black")) + theme(legend.position="none") +
  scale_y_continuous(name ="-log10 adjusted p-value")

dev.off()