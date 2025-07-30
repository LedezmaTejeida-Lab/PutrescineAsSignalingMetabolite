library(dplyr)
library(tidyr)
library(optparse)


option_list <- list(
    make_option(c("-t", "--TFname_PSSM_path"), type="character", default=NULL, 
              help="path of the TFname_PSSMv4.tsv file", metavar="character"),
    make_option(c("-i", "--TFinfo_path"), type="character", default=NULL, 
              help="path of the TF_info_v4.txt file", metavar="character"),
    make_option(c("-g", "--geneinfo_path"), type="character", default=NULL, 
              help="path of the geneannotation.tsv file", metavar="character"),
    make_option(c("-p", "--PSSM_info_path"), type="character", default=NULL, 
              help="path of thePSSM_info.tmp file", metavar="character"),
    make_option(c("-o", "--outpath"), type="character", default=NULL, 
              help="path of the PSSM_info.tsv file containing TFname, TFid, GeneCodifying, bnumber, No_sites, and width", metavar="character")
    
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$TFname_PSSM_path)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (TFname_PSSMv4.tsv file)", call.=FALSE)
}


# TFname_PSSM_path    <- "/space24/PGC/emhernan/masters/semester_2/PSSM_analysis_and_PhF/docs/TFname_PSSMv4.tsv"
# TFinfo_path         <- "/space24/PGC/emhernan/masters/semester_2/PSSM_analysis_and_PhF/docs/TF_info_v4.txt"
# geneinfo_path       <- "/space24/PGC/emhernan/masters/semester_2/PSSM_analysis_and_PhF/docs/geneannotation.tsv"
# PSSM_info_path      <- "/space24/PGC/emhernan/masters/semester_2/PSSM_analysis_and_PhF/docs/PSSM_info.tmp"
# outpath             <- "/space24/PGC/emhernan/masters/semester_2/PSSM_analysis_and_PhF/docs/PSSM_info.tsv"

TFname_PSSM <- read.table(file = opt$TFname_PSSM_path, sep = "\t", header = FALSE, col.names =c("TFname"))
TFinfo <- read.table(file = opt$TFinfo_path, sep = "\t", header = FALSE, col.names =c("TFid", "TFname", "GeneCodifying"))
geneinfo <- read.table(file = opt$geneinfo_path, sep = "\t", header = FALSE, col.names =c("geneName", "bnumber"), quote = "")
PSSM_info <- read.table(file = opt$PSSM_info_path, sep = "\t", header = FALSE, col.names =c("TFname", "No_sites", "width"))


df <- TFname_PSSM %>%
    left_join(TFinfo) %>%
    separate_rows(GeneCodifying, sep = ", ") %>%
    left_join(geneinfo, by = c("GeneCodifying"="geneName")) %>%
    group_by(TFname, TFid) %>%
    summarise(
    GeneCodifying = paste(GeneCodifying, collapse = ";"),
    bnumber = paste(bnumber, collapse = ";")
    ) %>%
    left_join(PSSM_info,by = c("TFname"= "TFname"))


write.table(x = df, file = opt$outpath, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# Four inconsistences:
# ECK120034977, ECK125257187, ECK120048948, ECK125257190 in PSSM file but other IDs in TFSet


# Transcription Factor ID: ECK120034977
# Transcription Factor Name: DcuR
# ECK125286588  DcuR    DNA-binding transcriptional activator DcuR

# Transcription Factor ID: ECK125257187
# Transcription Factor Name: Fur
# ECK125285350  Fur     DNA-binding transcriptional dual regulator Fur

# Transcription Factor ID: ECK120048948
# Transcription Factor Name: MetJ
# ECK125302602  MetJ    DNA-binding transcriptional repressor MetJ

# Transcription Factor ID: ECK125257190
# Transcription Factor Name: MntR
# ECK125302600  MntR    DNA-binding transcriptional dual regulator MntR,

