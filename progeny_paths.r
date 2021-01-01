#!/usr/bin/env Rscript

scriptDescription <- "A script that explores perturbation of canonical pathways via PROGENy."

scriptMandatoryArgs <- list(
  readCounts = list(
    abbr="-i",
    type="table",
    readoptions=list(sep="\t", stringsAsFactors=FALSE),
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
  clusterSamples = list(
    default=TRUE,
    type="logical",
    help="If samples should be clustered on the heatmap."
  ),
  expSpecies = list(
    default="Human",
    help="Species. Can be Human or Mouse."
  ),
  conditionColors = list(
    default=NULL,
    type="vector",
    help="Colors for points belonging to a given condition on PCA."
  ),
  commandRpath = list(
    default="commandR.r",
    help="Path to command line connectivity script (if not in cwd)."
  )
)

opt <- list()
for (rn in names(scriptOptionalArgs)){
  opt[[rn]] <- scriptOptionalArgs[[rn]][["default"]]
}

for (pk in c("tidyr", "dplyr", "purrr", "tibble", "DESeq2", "progeny", "pheatmap", "ggplot2", "cowplot", "ggplotify")){
  if(!(pk %in% (.packages()))){
    library(pk, character.only=TRUE)
  }
}

#' The main function of the script, executed only if called from command line.
#' Calls subfunctions according to supplied command line arguments.
#' 
#' @param opt list. a named list of all command line options; will be passed on 
#' 
#' @return Not intended to return anything, but rather to save outputs to files.
main <- function(opt){
  
  outFile <- opt$outFile
  opt$outFile <- NULL
  opt$help <- NULL
  opt$verbose <- NULL

  cat("Testing DE with EdgeR\n")
  progeny_results <- do.call(progenyPathwayScores, opt)
  cat("Saving tables\n")
  genetab2tsv(progeny_results$table, paste0(outFile, ".tsv"))
  
  cat("Saving figures\n")
  p <- plot_grid(progeny_results$figures$heat, progeny_results$figures$box, ncol=1)
  fig2pdf(p, paste0(outFile, ".pdf"), height=8.64, width=7.2)

  invisible(NULL)
}


#' Run PROGENy on a gene expression matrix to predict activity of canonical pathways.
#' 
#' @param readCounts data.frame. The count matrix; first column becoming the index, second labels.
#' @param conditionLabels character vector. Experimental condition labels. 
#' @param expSpecies string. Species, currently Human or Mouse. 
#' @param conditionOrder character vector. Order of conditions. Sensible to make control first. 
#' @param conditionColors character vector. Color associated to each experimental condition. 
#' @param clusterSamples If samples should be clustered on the heatmap.
#' @param ... Arguments inherited from command line but not used by this function. 
#' 
#' @return Pathway activity scores.
progenyPathwayScores <- function(readCounts, conditionLabels, expSpecies, conditionOrder=NULL, conditionColors=NULL, clusterSamples=TRUE, ...){

  ### Start by parsing inputs and setting some defaults
  if(is.null(conditionOrder)){
    conditionOrder = unique(conditionLabels)
  }
    
  conditions <- factor(conditionLabels, levels=conditionOrder)
  
  if (is.null(conditionColors)){
    conditionColors <- default_colors
  }

  metaData <- readCounts %>%
  colnames() %>%
  .[c(-1, -2)] %>%
  data.frame(group = conditions) %>%
  column_to_rownames(var = ".")

  ### PROGENy accepts a variance stabilized DEseq2 object as input
  gexInfo <- readCounts %>% 
    mutate(Gene = .[,2]) %>%
    .[,c(-1, -2)] %>%
    group_by(Gene) %>%
    summarise_all(mean) %>%
    ungroup() %>%
    column_to_rownames(var="Gene") %>%
    round() %>%
    DESeqDataSetFromMatrix(colData = metaData, design = ~group) %>%
    estimateSizeFactors() %>%
    estimateDispersions() %>%
    getVarianceStabilizedData() %>%
    progeny(scale=FALSE, organism=expSpecies) %>%
    data.frame() %>%
    rownames_to_column(var = "SampleID")
  
  ### Plot pathway perturbation scores with PROGENy

  matrixMax <- max(gexInfo[,-1])
  matrixMin <- min(gexInfo[,-1])
  matrixRng <- max(matrixMax, abs(matrixMax))
  matrixBrk <- pretty(-matrixRng:matrixRng, n = 10)

  p1 <- gexInfo %>%
  column_to_rownames(var = "SampleID") %>%
  t() %>%
  pheatmap(
      annotation_col=metaData,
      annotation_colors=list(group=setNames(conditionColors, conditionOrder)),
      annotation_legend=FALSE,
      annotation_names_row=FALSE,
      show_colnames=TRUE,
      show_rownames=TRUE,
      treeheight_col=0,
      cluster_cols=clusterSamples,
      color=colorRampPalette(c("#0000FF", "#1E90FF", "grey", "#FFA07A", "#FF0000"))(length(matrixBrk)),
      breaks=matrixBrk,
      display_numbers=TRUE,
      number_format="%.0f",
      main="PROGENy pathway activity scores"
      ) %>%
    as.ggplot()

  p2 <- gexInfo %>%
    pivot_longer(-SampleID, names_to="Pathway") %>%
    left_join(rownames_to_column(metaData, var="SampleID")) %>%
    ggplot(aes(x=Pathway, y=value, color=group, fill=group)) + 
    geom_boxplot(position=position_dodge2(width=0.8, preserve="single", padding=0.2), fill="white", outlier.shape=NA) +
    geom_jitter(position=position_jitterdodge(jitter.width=0.4, dodge.width=0.8), size=.5) +
    scale_color_manual(values=conditionColors, drop=FALSE) +
    theme_bw() +
    theme(
      axis.ticks = element_blank(),
      axis.text.x = element_text(size=10, angle=30, hjust=1),
      legend.position = "right",
      legend.title=element_blank()
    ) + 
    labs(x="", y="PROGENy score")

  invisible(list(figures=list(heat=p1, box=p2), table=gexInfo))

}

# Ensuring command line connectivity by sourcing an argument parser
source(opt$commandRpath, local=TRUE)