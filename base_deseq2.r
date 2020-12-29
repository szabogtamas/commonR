#!/usr/bin/env Rscript

scriptDescription <- "A script that runs a default DE anlysis with DESeq2 for Bulk RNAseq."

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
  labelVolcano = list(
    default=TRUE,
    type="logical",
    help="If the DE genes should be labelled on the Volcano plot."
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

for (pk in c("tidyr", "dplyr", "purrr", "tibble", "openxlsx", "DESeq2", "pheatmap", "EnhancedVolcano", "ggplot2", "cowplot", "ggplotify")){
  if(!(pk %in% (.packages()))){
    library(pk, character.only=TRUE)
  }
}

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

  cat("Testing DE with DESeq2\n")
  test_results <- do.call(testDEwithDeseq, opt)
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
  write.xlsx(test_results$tables, paste0(outFile, "summary.xlsx"))
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


#' Runs DE analysis with DESeq2 on a count matrix. Condition labels should be in the
#' order as samples appear in columns.
#' 
#' @param readCounts data.frame. The count matrix; first column becoming the index, second labels.
#' @param conditionLabels character vector. Experimental condition labels. 
#' @param conditionOrder character vector. Order of conditions. Sensible to make control first. 
#' @param conditionColors character vector. Color associated to each experimental condition. 
#' @param clusterSamples If samples should be clustered on the heatmap.
#' @param labelVolcano If the DE genes should be labelled on the Volcano plot. 
#' @param ... Arguments inherited from command line but not used by this function. 
#' 
#' @return Hitlists.
testDEwithDeseq <- function(readCounts, conditionLabels, conditionOrder=NULL, conditionColors=NULL, shrinkFc=FALSE, clusterSamples=TRUE, labelVolcano=TRUE, ...){

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
  data.frame(group = conditions) %>%
  column_to_rownames(var = ".")

  ### Run DESeq2
  de_test <- readCounts %>% 
    #mutate(Gene = .[,2]) %>%
    mutate(Gene = .[,1]) %>% #To harmonize with current EdgeR appproach, ENSEMBL ID is used instead of Symbol
    .[,c(-1, -2)] %>%
    #group_by(Gene) %>%
    #summarise_all(mean) %>%
    #ungroup() %>%
    column_to_rownames(var = "Gene") %>%
    DESeqDataSetFromMatrix(colData = metaData, design = ~group) %>%
    DESeq()
  
  normalized_counts <- log2(counts(de_test, normalized=TRUE) + 1)

  ### Format results before saving
  if(shrinkFc){
    de_res <- de_test %>%
      resultsNames() %>%
      .[-1] %>%
      map(~lfcShrink(de_test, .x))
  } else {
    de_res <- de_test %>%
      list() %>%
      rep(length(conditionOrder) - 1) %>%
      setNames(conditionOrder[seq(2, length(conditionOrder))]) %>%
      imap(~results(.x, contrast=c("group", .y, conditionOrder[1])))
  }
   
  de_tables <- de_res %>%
    map(data.frame)%>%
    map(~arrange(.x, padj))

  ### Compile a plot for every case (condition)
  plots <- list(
    heatm = imap(de_tables, draw_summary_heatmap, conditions, normalized_counts, conditionColors, geneDict, clusterSamples),
    volcano = imap(de_tables, draw_summary_volcano, conditions, geneDict, labelVolcano),
    mdp = map(de_res, draw_summary_mdplot)
  ) %>%
  pmap(draw_summary_panel)
  
  plots[["an_overview"]] <- draw_overview_panel(de_test, conditions, conditionColors, normalized_counts)

  de_tables[["normalized_matrix"]] <- as.data.frame(normalized_counts)

  if(all(geneDict == names(geneDict))){
    geneDict = NULL
  }

  invisible(list(figures=plots, tables=de_tables, rawResults=de_test, geneDict=geneDict))

}


#' Create an overview figure with a summary heatmap and a PCA plot.
#' 
#' @param de_test DESeqDataSet. Test result object from DESeq2, containing counts and metadata.
#' @param conditions named vector. All the conditions corresponding to matrix columns.
#' @param conditionColors named vector. Color codes for conditions. 
#' @param normalized_counts matrix. Counts normalized as log2(median of ratios); computed from de_test if not supplied.
#' 
#' @return overview plot woth PCA and heatmap.
draw_overview_panel <- function(de_test, conditions, conditionColors, normalized_counts=NULL){

  p1 <- de_test %>%
    DESeqTransform() %>%
    plotPCA(intgroup = "group") +
    scale_color_manual(values=conditionColors, drop=FALSE) +
    theme_bw()
  
  if(is.null(normalized_counts)){
    normalized_counts <- log2(counts(de_test, normalized=TRUE) + 1)
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
      #breaks=opt_color_breaks,
      color=colorRampPalette(c("white", "grey", "red"))(10)
      ) %>%
    as.ggplot()

  plot_grid(p1, p2, labels="AUTO", nrow=2, rel_heights=c(2, 1))
  
}


#' Create a multifigure panel with volcano and MD on the left and heatmap on the right.
#' 
#' @param heatm ggplot. Heatmap of top genes in the comparison.
#' @param volcano ggplot. Volcano plot summary of the comparison.
#' @param mdp ggplot. MD plot for QC purpose.
#' 
#' @return ggplotified DESeq2 MD plot.
draw_summary_panel <- function(heatm, volcano, mdp){
  
    plot_grid(
      plot_grid(volcano, mdp, labels="AUTO", nrow=2, rel_heights=c(2, 1)),
      heatm,
      ncol=2,
      labels=c("", "C"),
      hjust=0.2
    )

}


#' Draws an MD (mean vs. difference) plot for a given comparison.
#' 
#' @param de_test_result DESeq2 result. DE result of comparison for a single coef.
#' 
#' @return ggplotified DESeq2 MD plot.
draw_summary_mdplot <- function(de_test_result){

  de_test_result <<- de_test_result
  p <- expression(
    plotMA(de_test_result,
      main="",
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


#' Draws a summary volcano plot for showing change for all genes in a given comparison 
#' of a pair of conditions (one always being control).
#' 
#' @param result data.frame. Table of DE results after topTags for a single coef.
#' @param condition string. Experimental condition we want to compare to control. 
#' @param conditions factor. All the conditions, corresponding.
#' @param geneDict named vector. Names to be shown if renamed.
#' @param labelVolcano If the DE genes should be labelled on the Volcano plot. 
#' @param topNgene integer. How many genes should be shown on the Volcano plot.
#' @param ... arguments to be passed on to the EnhancedVolcano call. 
#' 
#' @return Hitlists.
draw_summary_volcano <- function(res, condition, conditions, geneDict=NULL, labelVolcano=TRUE, topNgene=25, ...){

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
  if(!labelVolcano){
    topResGen <- c()
    boxedLabels <- FALSE
    drawConnectors <- FALSE
  } else {
    boxedLabels <- TRUE
    drawConnectors <- TRUE
  }
  
  volcano.colors <- ifelse(
    res$padj < 0.05,
    ifelse(
      res$log2FoldChange < -2,
      "#3938fc",
      ifelse(res$log2FoldChange > 2, "#fa3316", "#808080")
    ),
    "#d3d3d3")
  names(volcano.colors) <- rownames(res)

  EnhancedVolcano(
    res,
    lab=res$gene,
    title=paste0('Differential Expression in\n', condition, ' vs. ', controllabel),
    titleLabSize = 12,
    x='log2FoldChange', y='padj',
    xlab = bquote(~Log2~ 'fold change'),
    pCutoff = 0.05,
    FCcutoff = 2,
    caption='', subtitle='',
    legendPosition='none',
    selectLab=topResGen,
    drawConnectors=drawConnectors, colConnectors="grey30", boxedLabels=boxedLabels,
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


#' Draws a summary heatmap of the top 50 genes for a pair of conditions (one always 
#' being control).
#' 
#' @param result data.frame. Table of DE results after topTags for a single coef.
#' @param condition string. Experimental condition we want to compare to control. 
#' @param conditions factor. All the conditions, corresponding columns in count matrix.
#' @param normalized_counts matrix. Normalized counts to be shown on the heatmap.
#' @param conditionColors named vector. Color codes for conditions. 
#' @param geneDict Mapping of human friendly gene names.
#' @param clusterSamples If samples should be clustered on the heatmap.
#' 
#' @return Hitlists.
draw_summary_heatmap <- function(res, condition, conditions, normalized_counts, conditionColors, geneDict, clusterSamples){
  
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
      cluster_cols = clusterSamples,
      color=colorRampPalette(c("white", "grey", "red"))(25),
      breaks=seq(0, 25, 1)
      #legend_labels = c(seq(0, 50, 5), "log2(Median of Ratios)\n")
      ) %>%
    as.ggplot()
}

# Ensuring command line connectivity by sourcing an argument parser
source(opt$commandRpath, local=TRUE)