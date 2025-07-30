#Name:
#	Confusion_matrix_PSSM.R
#Author:
#	Hernandez-Benitez Ericka M
#Version
#	v1
#Description
#	The script computes the components of a matrix confusion to evalute the predictive power of a PSSM
#	Input parameters:
#		--annotatedTFBS_path[file path]
#		--predictedTFBS_path[file path]
#		--output_file[file path]
#		--TF_name[character]
#   --FPR_TPR_table_path[file path]
#   --isStrongEvidence[boolean]
#

# annotatedTFBS_path  <- "/space24/PGC/emhernan/masters/semester_2/PSSM_analysis_and_PhF/docs/TF-TG-alltype-firstGene_locustag_cop.tsv"
# predictedTFBS_path  <- "/space24/PGC/emhernan/masters/semester_2/PSSM_analysis_and_PhF/rsat/matrix_scan/ArgP-predicted.tmp"
# output_file         <- "/space24/PGC/emhernan/masters/semester_2/PSSM_analysis_and_PhF/output/ArgP-metrics.tsv"
# TF_name             <- "ArgP"
# FPR_TPR_table_path  <- "/space24/PGC/emhernan/masters/semester_2/PSSM_analysis_and_PhF/output/FPR_TPR_table.tsv"
# isStrongEvidence <- FALSE
# -*- encoding: utf-8 -*-


print("...................................................Reading data...................................................")
# Loading libraries
library(dplyr)
library(stringr)
library(optparse)

option_list <- list(
    make_option(c("-a", "--annotatedTFBS_path"), type="character", default=NULL, 
              help="path of the TF-TG-firstOperonGene_locustag.tsv file", metavar="character"),
    make_option(c("-p", "--predictedTFBS_path"), type="character", default=NULL, 
              help="path of the predicted interactions (matrix_scan tmp file)", metavar="character"),
    make_option(c("-o", "--output_file"), type="character", default=NULL, 
              help="path of the output file where TF metrics should be saved", metavar="character"),
    make_option(c("-t", "--TF_name"), type="character", default=NULL, 
              help="TF name", metavar="character"),
    make_option(c("-m", "--FPR_TPR_table_path"), type="character", default=NULL, 
              help="path of the FPR_TPR_table.tsv file containing TFname, FPR and TPR", metavar="character"),
    make_option(c("-f", "--isStrongEvidence"), type="character", default=NULL, 
              help="boolean indicating if the evaluation must be performed with interactions with strong evidence", metavar="character")
    
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$annotatedTFBS_path)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (annotatedTFBS file path)", call.=FALSE)
}



# Loading data
annotatedTFBS <- read.table(file = opt$annotatedTFBS_path, header = TRUE, sep = "\t")
predictedTFBS <- read.table(file = opt$predictedTFBS_path, header = FALSE, sep = "\t")
isStrongEvidence <- opt$isStrongEvidence == "TRUE"

print("...................................................Done!...................................................")

print("...................................................Processing data...................................................")
# Processing data

colnames(predictedTFBS) <- c("TG_locustag", "TG_name", "ft_type")

predictedTFBS <- predictedTFBS %>%
                 group_by(TG_locustag, TG_name) %>%
                 filter(n() == 1 | ft_type == "site") %>%
                 ungroup()
                 
annotatedTFBS <- annotatedTFBS %>%
                  filter(tfName == opt$TF_name)

if(isStrongEvidence){
    annotatedTFBS <- annotatedTFBS %>% filter(riEvidence != "W") %>%
                        select(tfName, tgCopLocusTag) %>%
                        distinct()
    
} else{
    annotatedTFBS <- annotatedTFBS %>%
                        select(tfName, tgCopLocusTag) %>%
                        distinct()
}

TP <- annotatedTFBS %>%
        left_join(predictedTFBS, by = c("tgCopLocusTag" = "TG_locustag")) %>%
        filter(ft_type == "site") %>%
        nrow
FP <-  predictedTFBS %>%
        filter(TG_locustag %in% setdiff(predictedTFBS$TG_locustag, annotatedTFBS$tgCopLocusTag)) %>%
        filter(ft_type == "site") %>%
        nrow
FN <- annotatedTFBS %>%
        left_join(predictedTFBS, by = c("tgCopLocusTag" = "TG_locustag")) %>%
        filter(ft_type != "site" | is.na(ft_type)) %>%
        nrow
TN <- predictedTFBS %>%
        filter(TG_locustag %in% setdiff(predictedTFBS$TG_locustag, annotatedTFBS$tgCopLocusTag)) %>%
        filter(ft_type != "site") %>%
        nrow

confusion_matrix <- matrix(
                            c(TP, FP, FN, TN), 
                            nrow = 2, 
                            byrow = TRUE, 
                            dimnames = list(
                                c("Predicted (1)", "Predicted (0)"),
                                c("Actually (1)", "Actually (0)")
                                )
                            )


TPR <- TP/(TP +FN)
FPR <- FP/(FP+TN)

archive_name <- str_split(basename(opt$predictedTFBS_path), "-")[[1]][1]
parameters <- str_split(archive_name, "_")[[1]]

basename(opt$predictedTFBS_path)
print("...................................................Done!...................................................")

print("...................................................Saving data...................................................")
# Saving files
sink(opt$output_file)
    cat(sprintf('# TF name: %s with paramenters taxon_reference = %s, type_evidende = %s, min_w = %s , max_w = %s\n',
    opt$TF_name, parameters[1], parameters [2], parameters[3], parameters [4]))
    cat("# FPR\tTRP\n")
    cat("\tActually (1) \tActually (0)\n") 
sink()

write.table(x = confusion_matrix, file = opt$output_file , sep = "\t", col.names = FALSE, row.names = TRUE, quote = FALSE, append = TRUE)


sink(opt$FPR_TPR_table_path, append = TRUE)
    cat(paste(archive_name, FPR, TPR, sep = "\t"), "\n")
sink()


print("...................................................Done!...................................................")