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
  ),
  conditionColors = list(
    default=NULL,
    type="vector",
    help="Colors for points belonging to a given condition on PCA."
  )
)

opt <- list()
for (rn in names(scriptOptionalArgs)){
  opt[[rn]] <- scriptOptionalArgs[[rn]][["default"]]
}

for (pk in c("tidyr", "dplyr", "purrr", "tibble", "org.Hs.eg.db", "edgeR", "pheatmap", "ggplot2", "cowplot", "ggplotify")){
  if(!(pk %in% (.packages()))){
    library(pk, character.only=TRUE)
  }
}

default_colors <- c(
  '#1a476f', '#90353b', '#55752f', '#e37e00', '#6e8e84', '#c10534',
  '#938dd2', '#cac27e', '#a0522d', '#7b92a8', '#2d6d66', '#9c8847',
  '#bfa19c', '#ffd200', '#d9e6eb'
)

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
testDEwithEdgeR <- function(readCounts, conditionLabels, conditionOrder=NULL, conditionColors=NULL){
  
  #' Runs DE analysis with edger on a count matrix. Condition labels should be in the
  #' order as samples appear in columns.
  #' 
  #' @param readCounts data.frame. The count matrix; first column becoming the index.
  #' @param conditionLabels character vector. Experimental condition labels. 
  #' @param conditionOrder character vector. Order of conditions. Sensible to make control first. 
  #' 
  #' @return Hitlists.

  if(is.null(conditionOrder)){
    conditionOrder = unique(conditionLabels)
  }
  conditions <- factor(conditionLabels, levels=conditionOrder)
  design <- model.matrix(~conditions)

  if (is.null(conditionColors)){
    conditionColors <- default_colors
  }
  conditionColors <- conditionColors[
    rep(
      seq(1, length(conditionColors)),
      (length(conditionOrder)/length(conditionColors))+1
    )
  ]
  
  ### Use EdgeR for DE analysis
  y <- readCounts %>%
    column_to_rownames(colnames(.)[1]) %>%
    DGEList(group=conditions) %>%
    calcNormFactors() %>%
    estimateDisp(design)

  de_test <- y %>%
    glmFit(design) %>%
    list() %>%
    rep(length(conditionOrder) - 1) %>%
    map2(seq(2, length(conditionOrder)), glmLRT) %>%
    map(topTags, n="Inf") %>%
    map(pluck, 'table') %>%
    setNames(conditionOrder[seq(2, length(conditionOrder))])

}

# Ensuring command line connectivity by sourcing an argument parser
source("commandR.r")