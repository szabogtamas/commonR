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
  #' @param readCounts data.frame. The count matrix; first column becoming the index, second labels.
  #' @param conditionLabels character vector. Experimental condition labels. 
  #' @param conditionOrder character vector. Order of conditions. Sensible to make control first. 
  #' 
  #' @return Hitlists.

  if(is.null(conditionOrder)){
    conditionOrder = unique(conditionLabels)
  }
  conditionDict <- setNames(conditionLabels, colnames(readCounts[-1]))

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
  
  conditions <- factor(conditionLabels, levels=conditionOrder)
  design <- model.matrix(~conditions)
  geneDict <- setNames(readCounts[[2]], readCounts[[1]])

  ### Use EdgeR for DE analysis

  y <- readCounts %>%
      column_to_rownames(colnames(.)[1]) %>%
      .[,-1] %>%
      DGEList(group=conditions) %>%
      calcNormFactors() %>%
      estimateDisp(design)
  
  normalized_counts <- y %>%
  cpm()+1 %>%
  log2()

  de_test <- y %>%
    glmFit(design) %>%
    list() %>%
    rep(length(conditionOrder) - 1) %>%
    map2(seq(2, length(conditionOrder)), glmLRT) %>%
    #map(topTags, n="Inf") %>%
    #map(pluck, 'table') %>%
    setNames(conditionOrder[seq(2, length(conditionOrder))])
  
  heatmaps <- map2(de_test, names(de_test), draw_summary_heatmap, conditions, normalized_counts, conditionColors)
  volcanos <- map2(de_test, names(de_test), draw_summary_volcano, conditions, geneDict)
  mdplots <- map(
    de_test,
    plotMD, 
    main="",
    legend=FALSE,
    bty="L", 
    ex=0.8,
    cex.lab=0.8,
    cex.axis=0.8,
    mgp=c(1.5,0.5,0)
  ) %>%
    map(expression) %>%
    map(
      as.ggplot(vjust=0, scale=1) +
      theme(plot.margin = unit(c(0,0,0,-0.1), "in"))
    )

  plots <- map(seq(1, length(conditionOrder)) %>%
    plot_grid(plot_grid(p1, p2, p3, labels="AUTO", nrow=3, rel_heights=c(2, 1, 1)), p4, ncol=2, labels=c("", "D"), hjust=0.2)

  p1 <- pheatmap(normalized_counts)
  p2 <- as.ggplot(expression(plotMDS(y, col=conditionColors[conditions])))
  
  plots[["Summary"]] <- plot_grid(p1, p2, labels="AUTO", nrow=2, rel_heights=c(2, 1))
}


draw_summary_volcano <- function(de_test_result, condition, conditions, geneDict=NULL, topNgene=25, ...){

  #' Draws a summary volcano plot for showing change for all genes in a given comparison 
  #' of a pair of conditions (one always being control).
  #' 
  #' @param de_test_result edgeR result. DE result of comparison for a single coef.
  #' @param condition string. Experimental condition we want to compare to control. 
  #' @param conditions factor. All the conditions, corresponding.
  #' @param geneDict named vector. Names to be shown if renamed.
  #' @param topNgene integer. How many genes should be shown on the Volcano plot. 
  #' @param ... arguments to be passed on to the EnhancedVolcano call. 
  #' 
  #' @return Hitlists.
  
  res <- topTags(de_test_result, n="Inf")$table
  if(is.null(geneDict)){
    res$gene <- rownames(res)
  } else {
    res$gene <- as.character(geneDict[rownames(res)])
  }
  topResGen <- res$gene %>%
    .[.!=""] %>%
    unique() %>%
    .[1:topNgene]
  controllabel <- levels(conditions)[[1]]
  
  volcano.colors <- ifelse(
    res$FDR < 0.05,
    ifelse(
      res$logFC < -2,
      "#3938fc",
      ifelse(res$logFC > 2, "#fa3316", "#808080")
    ),
    "#d3d3d3")
  names(volcano.colors) <- volcano.colors
  EnhancedVolcano(
    res,
    lab=res$gene,
    title=paste0('Differential Expression in\n', condition, ' vs. ', controllabel),
    titleLabSize = 12,
    x='logFC', y='FDR',
    xlab = bquote(~Log~ 'fold change'),
    pCutoff = 0.05,
    FCcutoff = 2,
    caption='', subtitle='',
    legendPosition='none',
    selectLab=topResGen,
    drawConnectors=TRUE, colConnectors="grey30", boxedLabels=TRUE,
    gridlines.major=FALSE, gridlines.minor=FALSE,
    #colCustom = volcano.colors,
    borderWidth = 0.3,
    axisLabSize=10, subtitleLabSize=1, captionLabSize=1,
    ...
    ) + 
  theme(
    axis.text.x = element_text(color="black", size=10),
    axis.text.y = element_text(color="black", size=10),
    plot.margin = unit(c(0.1,0,-0.1,0.25), "in")
    )
}


draw_summary_heatmap <- function(de_test_result, condition, conditions, normalized_counts, conditionColors, geneDict){
  
  #' Draws a summary heatmap of the top 50 genes for a pair of conditions (one always 
  #' being control).
  #' 
  #' @param de_test_result edgeR result. DE result of comparison for a single coef.
  #' @param condition string. Experimental condition we want to compare to control. 
  #' @param conditions factor. All the conditions, corresponding columns in count matrix.
  #' @param normalized_counts matrix. Normalized counts to be shown on the heatmap.
  #' @param conditionColors named vector. Color codes for conditions. 
  #' 
  #' @return Hitlists.
  
  res <- topTags(de_test_result, n="Inf")
  samples_to_show <- which(conditions == condition | conditions == levels(conditions)[[1]])
  annots <- data.frame(condition=as.character(conditions[samples_to_show]))
  rownames(annots) <- colnames(normalized_counts)[samples_to_show]
  
  normalized_counts  %>%
    .[rownames(res[1:50,]),samples_to_show] %>%
    pheatmap(
      annotation_col=annots,
      annotation_colors=list(condition=conditionColors),
      annotation_legend=FALSE,
      labels_row=geneDict[rownames(res[1:50,])],
      show_colnames=TRUE,
      cluster_rows=FALSE,
      color=colorRampPalette(c("white", "grey", "red"))(10),
      breaks=seq(0, 50, 5)
      ) %>%
    as.ggplot()
}
}


### Patching Pheatmap with diagonal column labels;
### Idea from https://www.thetopsites.net/article/54919955.shtml
draw_colnames_30 <- function (coln, gaps, ...) {
    coord = pheatmap:::find_coordinates(length(coln), gaps)
    x = coord$coord - 0.5 * coord$size
    res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 30, gp = gpar(...))
    return(res)}

assignInNamespace(x="draw_colnames", value="draw_colnames_30", ns=asNamespace("pheatmap"))

# Ensuring command line connectivity by sourcing an argument parser
source("commandR.r")