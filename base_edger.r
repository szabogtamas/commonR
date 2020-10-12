#!/usr/bin/env Rscript

scriptDescription <- "A script that runs a default DE anlysis with edgeR for Bulk RNAseq."

scriptMandatoryArgs <- list(
  readCounts = list(
    abbr="-i",
    type="table",
    readoptions=list(sep="\t", stringsAsFactors=FALSE)
    help="Table of raw counts. A bit cleaned up: first column is GeneIDs, rest is counts."
  ),
  conditionLabels = list(
    abbr="-l",
    type="vector",
    help="Labels for samples (same order as count matrix columns!), reflecting experimental conditions."
  ),
  outFile = list(
    abbr="-o",
    help="Prefix for output files."
  )
)

scriptOptionalArgs <- list(
  conditionOrder = list(
    default=NULL,
    type="vector",
    help="Order of conditions in the experimental design formula. Makes sense to put control as first."
  )
)

opt <- list()
for (rn in names(scriptOptionalArgs)){
  opt[[rn]] <- scriptOptionalArgs[[rn]][["default"]]
}

for (pk in c("tidyr", "dplyr", "edgeR", "ggplot2", "cowplot", "ggplotify")){
  if(!(pk %in% (.packages()))){
    library(pk, character.only=TRUE)
  }
}


### Read command line parameters or use defaults

countfile <- "/home/szabo/myScratch/CrisprScreen/pipelinetesting/intermediates/15_bam_sequences/read_counts_guidelevel.txt"

controlbars <- "GCATGC,CTAGCA,CGTACG,ACTGAC"
casebars <- "TAGCTA,GACTGA,GATCGT,CTGACT"

controllabel <- "Control"
caselabel <- "Case"

condition_colors <- "black,red"

### Define a main function that will only be executed if called from command line
main <- function(opt){
  
  #' The main function of the script, executed only if called from command line.
  #' Calls subfunctions according to supplied command line arguments.
  #' 
  #' @param opt list. a named list of all command line options; will be passed on 
  #' 
  #' @return Not intended to return enything, but rather save outputs to files.
  
  outFile <- opt$outFile
  opt$outFile <- NULL
  opt$help <- NULL

  testDE(opt)
  cat("Saving figure\n")
  pdf(paste0(outFile, ".pdf"), height=9.6, width=7.2)
  print(p)
  dev.off()

}

### The actual working horse, called by main()
testDE <- function(readCounts, conditionLabels, conditionOrder){
  
  #' The main function of the script, executed only if called from command line.
  #' Calls subfunctions according to supplied command line arguments.
  #' 
  #' @param opt list. a named list of all command line options; will be passed on 
  #' 
  #' @return Not intended to return enything, but rather save outputs to files.

  if(is.null(conditionOrder)){
    conditionOrder = unique(conditionLabels)
  }
  conditions <- factor(conditionLabels, levels=conditionOrder)
  design <- model.matrix(~conditions)

  ### Deal with colors later!!!
  condition_colors <- unlist(strsplit(condition_colors, ",", fixed=TRUE))
  names(condition_colors) <- levels(conditions)

  ### Use EdgeR for DE analysis
  y <- readCounts %>%
    column_to_rownames(1)
    DGEList(group=conditions) %>%
    calcNormFactors() %>%
    estimateDisp(design)

  de_test <- y %>%
    glmFit(design) %>%
    glmLRT()

  ### Make this a more general case, returning hitlists for all conditions
  res <- topTags(lrt, n="Inf")$table

}

# Ensuring command line connectivity by sourcing an argument parser
source("commandR.r")