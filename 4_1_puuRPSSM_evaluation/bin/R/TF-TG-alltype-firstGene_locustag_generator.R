# Loading libraries
library(dplyr)
library(optparse)

option_list <- list(
    make_option(c("-i", "--interactionPath"), type="character", default=NULL, 
              help="path of the TF-TG-firstOperonGene.tsv file", metavar="character"),
    make_option(c("-a", "--annotationPath"), type="character", default=NULL, 
              help="path of the gene-bnumber.tsv file", metavar="character"),
    make_option(c("-o", "--outPath"), type="character", default=NULL, 
              help="path of the TF-TG-firstOperonGene_locustag.tsv file", metavar="character")
    
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$interactionPath)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (interactionPath file)", call.=FALSE)
}


# interactionPath <- "/home/ericka/docs/TF-TG-alltype-firstGene.tsv"
# annotationPath <- "/home/ericka/docs/gene-bnumber.tsv"
# outPath <- "/home/ericka/docs/TF-TG-alltype-firstGene_locustag.tsv"

print(".................................. Loading data ..................................")
interactions <- read.table(opt$interactionPath, sep = "\t", header = FALSE)
annotation <- read.table(opt$annotationPath, sep = "\t", header = TRUE, quote = "")
print(".................................. Done! ..................................")


print(".................................. Processing data ..................................")
colnames(interactions) <- c("riType", "tfName","tgName", "riEvidence")

df <- interactions %>%
    left_join(annotation, by  = c("tgName"="geneName")) %>%
    rename(tgLocusTag = bnumber)

print(".................................. Done! ..................................")

print(".................................. Saving data ..................................")
write.table(df, opt$outPath, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
print(".................................. Done! ..................................")