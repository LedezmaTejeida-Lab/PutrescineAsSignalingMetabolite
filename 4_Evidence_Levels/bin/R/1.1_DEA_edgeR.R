## RNA-seq analysis with edgeR
## The script runs under the DEA R project

## Loading libraries
library(edgeR) 
library(dplyr)
library(tidyr)
library(tibble)
library(statmod)
library(RColorBrewer)


## Reading raw count tables
raw_counts <- read.table(file = "/home/emhernan/4_Evidence_Levels/input/Cuentas_RNAseq_Ch24-10-v2.tsv/Cuentas_RNAseq_Ch24-10-v2.tsv", header = TRUE, sep = "\t")

## Processing file
raw_counts_p <- raw_counts %>%
  select(Geneid, C.1, C.3, C.4, F.1, F.2, F.4) %>%
  rename(B.1 = F.1, B.2 = F.2, B.4 = F.4) %>%
  column_to_rownames(var = "Geneid") 

## Checking depth sequencing for each replicate (million reads per sample)       
head(raw_counts_p)
dim(raw_counts_p)
round(colSums(raw_counts_p)/1e6)

# C.1 C.3 C.4 F.1 F.2 F.4 M.1 M.3 M.4 H.1 H.2 H.3 
# 14  14  15  15  16  15  15  15  14  12  16  13 

## Normalize the raw count to count per million (CPM). 
## Filtering the gene that have low expression values in the majority of samples.
## A generally recommended cutoff of read number for a low-expressed transcript is CPM of 1
## In the case of a typical sequencing depth of a total 10–30 million reads per sample, this cutoff 
## corresponds to 10–30 reads mapped to the transcript.
## For example, only keep the genes whose CPM value is higher than 1 in at three samples

keep <- rowSums(cpm(raw_counts_p) >= 1) >= 3
table(keep) # 9 genes does not pass the cut-off

counts <- raw_counts_p[keep,]
dim(counts)

## Extracting group name
group <- factor(sub("..$", "", colnames(counts)))
group
table(group) ## double checking the number of replicates per condition



## Creating a Digital Gene Expression (DGEList object)
cds <- DGEList(counts=counts, group=group)
cds <- cds[rowSums(1e+06 * (cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)))> 1) >= 3, ] ## Double filtering of CPM > 1 in at least 3 conditions
keep2 <- filterByExpr(cds)
cds <- cds[keep2, ,keep.lib.sizes=FALSE] ## another filter for low expressed genes
cds <- calcNormFactors(cds, method = "TMM")

cds$samples$norm.factors
  
  
## Plotting MSD to check how replicates cluster

colors <- brewer.pal(n = 4, name = "Dark2")
names(colors) <- levels(group)
png(filename = "/home/emhernan/4_Evidence_Levels/png/MDSplot.png", width=7*300, height=7*300, res=300)
plotMDS(cds, col=colors[cds$samples$group], labels=colnames(cds$counts))
dev.off()


## Coding experimental design

design <- model.matrix(~0+cds$samples$group)
colnames(design) <- levels(cds$samples$group)
design

### Estimating dispersion
cds <- estimateDisp(cds, design=design, robust = TRUE)
plotBCV(cds)

fit <- glmQLFit(cds, design=design, robust=TRUE)
plotQLDisp(fit)

## Performing DEA contrasts

contrast <- makeContrasts(
  "B" = "B - C",
  levels=cds$design
)
contrast

qlf <- glmQLFTest(fit, contrast=contrast)

## Adjustment of p-values using the Benjamini-Hochberg method 
## Add the adjusted p-values to the results table
p_values <- qlf$table$PValue
qlf$table$adj.PVal <- p.adjust(p_values, method="BH")
topTable <- topTags(qlf, n=Inf)$table
head(topTable)



DEA_gene_results <- topTable %>% rownames_to_column(var = "Geneid") %>%
  left_join(raw_counts %>% select(Geneid, Locus, ID)) %>%
  filter(abs(logFC) > 1) %>%
  filter(FDR < 0.05 & adj.PVal < 0.05)

all_genes <- topTable %>% 
              rownames_to_column(var = "Geneid") %>%
              left_join(raw_counts %>% select(Geneid, Locus, ID))
write.table(x = DEA_gene_results, file = "/home/emhernan/4_Evidence_Levels/output/DEA_gene_results_bean.tsv",quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)
write.table(x = all_genes, file = "/home/emhernan/4_Evidence_Levels/output/all_gene_results_bean.tsv",quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)
