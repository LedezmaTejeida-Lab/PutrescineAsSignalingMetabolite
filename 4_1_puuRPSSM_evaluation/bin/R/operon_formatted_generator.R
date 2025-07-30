library(dplyr)
library(tidyr)
library(stringr)
library(optparse)


option_list <- list(
    make_option(c("-b", "--gene_bnumber_path"), type="character", default=NULL, 
              help="path of the gene-bnumber.tsv file", metavar="character"),
    make_option(c("-p", "--operons_path"), type="character", default=NULL, 
              help="path of the operon_genes.tsv file", metavar="character"),
    make_option(c("-o", "--outpath"), type="character", default=NULL, 
              help="path of the eco.ope file", metavar="character")
    
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$gene_bnumber_path)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (gene-bnumber.tsv file)", call.=FALSE)
}

# gene_bnumber_path   <- "C:/Users/monts/Desktop/gene-bnumber.tsv"
# operons_path        <- "C:/Users/monts/Desktop/operon_genes.tsv"
# outpath             <- "C:/Users/monts/Desktop/eco.ope"


gene_bnumber <- read.table(file = opt$gene_bnumber_path, header = TRUE, sep = "\t", quote = "")
operons     <- read.table(file = opt$operons_path, header = FALSE, sep = "\t", quote = "", col.names = c("geneName"))

new_rows  <- tibble(
          geneName  = c("allF", "allG", "allH", "allK", "mogA", "eflP", "madA") , 
          bnumber = c("b0518", "b4572", "b0520", "b0521", "b0009", "b2171", "b3888"))

gene_bnumber <- gene_bnumber %>%
  add_row(new_rows)


df <- operons %>%
  mutate(operonID = 1:nrow(operons)) %>%
  separate_rows(geneName, sep = ";") %>%
  left_join(gene_bnumber) %>%
 mutate(bnumber = paste("eco-", bnumber, sep = "")) %>%
  group_by(operonID) %>%
  summarise(
    geneNames = str_c(geneName, collapse = " "),
    bnumbers = str_c(bnumber, collapse = " ")
  ) %>%
  ungroup() %>%
  select(bnumbers)

write.table(x = df, file = opt$outpath, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


# geneName operonID bnumber
# allF          870 b0518     
# allG          870 b4572     
# allH          870 b0520     
# allK          870 b0521     
# mogA         1004 b0009     
# eflP         1342 b2171     
# madA         2137 b3888
# sroG           37 None   
# insEF3         68 None   
# insCD6        114 None   
# insCD1        149 None   
# insAB3        167 None   
# tpke70        168 None   
# insAB6        221 None   
# insCD4        337 None   
# insAB1        471 None   
# sroD          581 None   
# insCD3        786 None   
# insEF5        813 None   
# sraA          849 None   
# insEF2       1177 None   
# insEF1       1240 None   
# insEF4       1314 None   
# insAB4       1364 None   
# gndA         1686 None   
# insCD2       1840 None   
# insCD5       1957 None   
# insAB5       1996 None   
# insAB2       2080 None   
# och5         2528 None   
# tpke11       2606 None 

# None: are trasnposes or IS that does not have a locus-tag
# NA:synonims problem
