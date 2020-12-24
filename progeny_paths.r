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

for (pk in c("tidyr", "dplyr", "purrr", "tibble", "openxlsx", "edgeR", "pheatmap", "EnhancedVolcano", "ggplot2", "cowplot", "ggplotify")){
  if(!(pk %in% (.packages()))){
    library(pk, character.only=TRUE)
  }
}
myScratch/DiffExpr/rnaSeqEmuMyc/tables/counts.tsv
#TODO This script will use Progeny to show which pathways are perturbed. Chack out https://github.com/saezlab/progeny/blob/master/vignettes/progenyBulk.Rmd

#' The main function of the script, executed only if called from command line.
#' Calls subfunctions according to supplied command line arguments.
#' 
#' @param opt list. a named list of all command line options; will be passed on 
#' 
#' @return Not intended to return anything, but rather save outputs to files.
main <- function(opt){
  
  outFile <- opt$outFile
  opt$outFile <- NULL
  opt$help <- NULL
  opt$verbose <- NULL

  cat("Testing DE with EdgeR\n")
  progeny_results <- do.call(progenyPathwayScores, opt)
  cat("Saving tables\n")
  geneDict <- progeny_results$geneDict
  genetab2tsv(progeny_results$table, outFile)
  
  cat("Saving figures\n")
  fignames <- test_results$figures %>%
    names() %>%
    gsub("/", "___", .) %>%
    paste0(outFile, .)
  map2(
    test_results$figures,
    fignames,
    fig2pdf,
    height=8.64,
    width=7.2
  )

  invisible(NULL)
}


#' Runs DE analysis with edger on a count matrix. Condition labels should be in the
#' order as samples appear in columns.
#' 
#' @param readCounts data.frame. The count matrix; first column becoming the index, second labels.
#' @param conditionLabels character vector. Experimental condition labels. 
#' @param expSpecies string. Species, currently Human or Mouse. 
#' @param conditionOrder character vector. Order of conditions. Sensible to make control first. 
#' @param conditionColors character vector. Color associated to each experimental condition. 
#' @param clusterSamples If samples should be clustered on the heatmap.
#' @param ... Arguments inherited from command line but not used by this function. 
#' 
#' @return Hitlists.
progenyPathwayScores <- function(readCounts, conditionLabels, expSpecies, conditionOrder=NULL, conditionColors=NULL, clusterSamples=TRUE, ...){

  ### Start by parsing inputs and setting some defaults
  if(is.null(conditionOrder)){
    conditionOrder = unique(conditionLabels)
  }
  conditionDict <- setNames(conditionLabels, colnames(readCounts[c(-1, -2)]))
    
  conditions <- factor(conditionLabels, levels=conditionOrder)
  
  if (is.null(conditionColors)){
    conditionColors <- default_colors
  }
  conditionColors <- conditionColors[
    rep(
      seq(1, length(conditionColors)),
      (length(conditionOrder)/length(conditionColors))+1
    )
  ]
  conditionColors <- setNames(conditionColors[1:length(levels(conditions))], levels(conditions))
  
  geneDict <- setNames(readCounts[[2]], readCounts[[1]])

  metaData <- readCounts %>%
  colnames() %>%
  .[c(-1, -2)] %>%
  data.frame(group = conditionLabels) %>%
  column_to_rownames(var = ".")

  ### PROGENy accepts a variance stabilized DEseq2 object as input
  gexInfo <- readCounts %>% 
    mutate(Gene = .[,2]) %>%
    .[,c(-1, -2)] %>%
    group_by(Gene) %>%
    summarise_all(mean) %>%
    ungroup() %>%
    column_to_rownames(var = "Gene") %>%
    DESeqDataSetFromMatrix(colData = metaData, design = ~group) %>%
    estimateSizeFactors() %>%
    estimateDispersions() %>%
    getVarianceStabilizedData() %>%
    progeny(scale=FALSE, organism = expSpecies) %>%
    data.frame() %>%
    rownames_to_column(var = "SampleID")
  
  ### Plot pathway perturbation scores with PROGENy
  p1 <- gexInfo %>%
    pivot_longer(-SampleID, names_to = "Pathway") %>%
    left_join(rownames_to_column(metaData, var = "SampleID")) %>%
    ggplot(aes(x=Pathway, y=value, color=group, fill=group)) + 
    geom_boxplot(position=position_dodge2(width = 0.8, preserve = "single", padding = 0.2), fill = "white", outlier.shape = NA) +
    geom_jitter(position=position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8)) +
    scale_color_manual(values=conditionColors, drop=FALSE) +
    theme(
      axis.ticks = element_blank(),
      axis.text.x = element_text(size=7, angle=30, hjust=1),
      legend.position = "right"
    ) + 
    labs(x="", y="PROGENy score")

  p2 <- gexInfo %>%
    pivot_longer(-SampleID, names_to = "Pathway") %>%
    left_join(rownames_to_column(metaData, var = "SampleID")) %>%
    ggplot(aes(x=Pathway, y=value, color=group, fill=group)) + 
    geom_boxplot(position=position_dodge2(width = 0.8, preserve = "single", padding = 0.2), fill = "white", outlier.shape = NA) +
    geom_jitter(position=position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8)) +
    scale_color_manual(values=conditionColors, drop=FALSE) +
    theme(
      axis.ticks = element_blank(),
      axis.text.x = element_text(size=7, angle=30, hjust=1),
      legend.position = "right"
    ) + 
    labs(x="", y="PROGENy score")

  p2 <- gexInfo %>%
  column_to_rownames(var = "SampleID") %>%
  t() %>%
  pheatmap(
      annotation_colors=list(condition=conditionColors),
      annotation_legend=FALSE,
      annotation_names_row=FALSE,
      show_colnames=FALSE,
      show_rownames=TRUE,
      treeheight_col=0,
      color=colorRampPalette(c("white", "grey", "red"))(10)
      ) %>%
    as.ggplot()

  invisible(list(figures=plots, table=gexInfo))

  annots <- data.frame(condition=as.character(conditions[colnames(normalized_counts)]))
  rownames(annots) <- colnames(normalized_counts)
  
  max_variances <- normalized_counts %>%
    apply(1, var) %>%
    order(decreasing=TRUE) %>%
    .[1:500]
  
  opt_color_breaks <- normalized_counts %>%
    .[max_variances,]  %>%
    quantile(seq(0.1, 1, 0.1)) %>%
    as.numeric()
    
  p2 <- normalized_counts %>%
  .[max_variances,]  %>%
  t() %>%
  pheatmap(
      annotation_row=annots,
      annotation_colors=list(condition=conditionColors),
      annotation_legend=FALSE,
      annotation_names_row=FALSE,
      show_colnames=FALSE,
      show_rownames=TRUE,
      treeheight_col=0,
      color=colorRampPalette(c("white", "grey", "red"))(10),
      breaks=opt_color_breaks
      ) %>%
    as.ggplot()

  normalized_counts  %>%
    .[rownames(res[1:50,]),samples_to_show] %>%
    pheatmap(
      annotation_col=annots,
      annotation_colors=list(condition=conditionColors),
      annotation_legend=FALSE,
      labels_row=geneDict[rownames(res[1:50,])],
      show_colnames=TRUE,
      cluster_rows=FALSE,
      cluster_cols = clusterSamples,
      color=colorRampPalette(c("white", "grey", "red"))(10),
      breaks=seq(0, 50, 5)
      ) %>%
    as.ggplot()
}

# Ensuring command line connectivity by sourcing an argument parser
source(opt$commandRpath, local=TRUE)