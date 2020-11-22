#!/usr/bin/env Rscript

scriptDescription <- "A script that runs a default DE anlysis with edgeR for Bulk RNAseq."

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

for (pk in c("tidyr", "dplyr", "purrr", "tibble", "edgeR", "pheatmap", "EnhancedVolcano", "ggplot2", "cowplot", "ggplotify")){
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
  opt$verbose <- NULL

  cat("Testing DE with EdgeR\n")
  test_results <- do.call(testDEwithEdgeR, opt)
  cat("Saving tables\n")
  geneDict <- test_results$geneDict
  tablenames <- test_results$tables %>%
    names() %>%
    gsub("/", "___", .) %>%
    paste0(outFile, .)
  map2(
    test_results$tables,
    tablenames,
    genetab2tsv,
    relabels=geneDict
   ) 
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

### The actual working horse, called by main()
testDEwithEdgeR <- function(readCounts, conditionLabels, conditionOrder=NULL, conditionColors=NULL, ...){
  
  #' Runs DE analysis with edger on a count matrix. Condition labels should be in the
  #' order as samples appear in columns.
  #' 
  #' @param readCounts data.frame. The count matrix; first column becoming the index, second labels.
  #' @param conditionLabels character vector. Experimental condition labels. 
  #' @param conditionOrder character vector. Order of conditions. Sensible to make control first. 
  #' @param conditionColors character vector. Color associated to each experimental condition. 
  #' @param ... Arguments inherited from command line but not used by this function. 
  #' 
  #' @return Hitlists.

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
    setNames(conditionOrder[seq(2, length(conditionOrder))])
  
  de_tables <- de_test %>% 
    map(topTags, n="Inf") %>%
    map(pluck, 'table')

  ### Compile a plot for every case (condition)
  plots <- list(
    heatm = map2(de_tables, names(de_tables), draw_summary_heatmap, conditions, normalized_counts, conditionColors, geneDict),
    volcano = map2(de_tables, names(de_tables), draw_summary_volcano, conditions, geneDict),
    mdp = map(de_test, draw_summary_mdplot)
  ) %>%
  pmap(draw_summary_panel)

  plots[["an_overview"]] <- draw_overview_panel(y, conditionDict, conditionColors, normalized_counts)

  de_tables[["normalized_matrix"]] <- as.data.frame(normalized_counts)

  if(all(geneDict == names(geneDict))){
    geneDict = NULL
  }

  invisible(list(figures=plots, tables=de_tables, geneDict=geneDict))

}


draw_overview_panel <- function(y, conditions, conditionColors, normalized_counts = NULL){

  #' Create an overview figure with a summary heatmap and a PCA plot.
  #' 
  #' @param y edgeR matrix. Count matrix in EdgeR after normalization and dispersion calculation.
  #' @param conditions named vector. All the conditions corresponding to matrix columns.
  #' @param conditionColors named vector. Color codes for conditions. 
  #' @param normalized_counts matrix. Counts normalized as logCPM; computed from y if not supplied.
  #' 
  #' @return ggplotified EdgeR MD plot.

  y <<- y
  mds_cls <<- conditionColors[conditions]
  mds_lvls <<- unique(names(conditionColors))
  p1 <- expression(
    plotMDS(y,
      col=mds_cls,
      main="",
      bty="L", 
      cex=0.8,
      cex.lab=0.8,
      cex.axis=0.8,
      mgp=c(1.5,0.5,0)
    ),
    legend(
      "bottomright",
      legend=mds_lvls, 
      col=mds_cls, 
      pch=19, 
      bty="n", 
      pt.cex=1, 
      cex=0.8,
      y.intersp=2
    )
  ) %>%
  as.ggplot(vjust=0, scale=1) +
    theme(plot.margin = unit(c(0,0,0,-0.1), "in"))
  
  if(TRUE){ #if(is.null(normalized_counts)){
    normalized_counts <- y %>%
      cpm()+1 %>%
      log2()
  }
  
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

  plot_grid(p1, p2, labels="AUTO", nrow=2, rel_heights=c(2, 1))
  
}


draw_summary_panel <- function(heatm, volcano, mdp){

  #' Create a multifigure panel with volcano and MD on the left and heatmap on the right.
  #' 
  #' @param heatm ggplot. Heatmap of top genes in the comparison.
  #' @param volcano ggplot. Volcano plot summary of the comparison.
  #' @param mdp ggplot. MD plot for QC purpose.
  #' 
  #' @return ggplotified EdgeR MD plot.
  
    plot_grid(
      plot_grid(volcano, mdp, labels="AUTO", nrow=2, rel_heights=c(2, 1)),
      heatm,
      ncol=2,
      labels=c("", "C"),
      hjust=0.2
    )

}


draw_summary_mdplot <- function(de_test_result){

  #' Draws an MD (mean vs. difference) plot for a given comparison.
  #' 
  #' @param de_test_result edgeR result. DE result of comparison for a single coef.
  #' 
  #' @return ggplotified EdgeR MD plot.
  
  de_test_result <<- de_test_result
  p <- expression(
    plotMD(de_test_result,
      main="",
      legend=FALSE,
      bty="L", 
      cex=0.8,
      cex.lab=0.8,
      cex.axis=0.8,
      mgp=c(1.5,0.5,0)
    )
  )
  as.ggplot(p, vjust=0, scale=1) +
    theme(plot.margin = unit(c(0,0,0,-0.1), "in"))
}


draw_summary_volcano <- function(res, condition, conditions, geneDict=NULL, topNgene=25, ...){

  #' Draws a summary volcano plot for showing change for all genes in a given comparison 
  #' of a pair of conditions (one always being control).
  #' 
  #' @param result data.frame. Table of DE results after topTags for a single coef.
  #' @param condition string. Experimental condition we want to compare to control. 
  #' @param conditions factor. All the conditions, corresponding.
  #' @param geneDict named vector. Names to be shown if renamed.
  #' @param topNgene integer. How many genes should be shown on the Volcano plot. 
  #' @param ... arguments to be passed on to the EnhancedVolcano call. 
  #' 
  #' @return Hitlists.
  
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
  names(volcano.colors) <- rownames(res)
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
    colCustom = volcano.colors,
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


draw_summary_heatmap <- function(res, condition, conditions, normalized_counts, conditionColors, geneDict){
  
  #' Draws a summary heatmap of the top 50 genes for a pair of conditions (one always 
  #' being control).
  #' 
  #' @param result data.frame. Table of DE results after topTags for a single coef.
  #' @param condition string. Experimental condition we want to compare to control. 
  #' @param conditions factor. All the conditions, corresponding columns in count matrix.
  #' @param normalized_counts matrix. Normalized counts to be shown on the heatmap.
  #' @param conditionColors named vector. Color codes for conditions. 
  #' 
  #' @return Hitlists.
  
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

# Ensuring command line connectivity by sourcing an argument parser
source(opt$commandRpath, local=TRUE)