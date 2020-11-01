#!/usr/bin/env Rscript

scriptDescription <- "A script that takes a DE result with numerical changes under a condition and runs GSEA."

scriptMandatoryArgs <- list(
  scoreTables = list(
    abbr="-i",
    type="tables",
    help="Tables with unique ID and name as first two column and numeric as rest."
  ),
  outFile = list(
    abbr="-o",
    help="Base name for output files."
  ),
  outPrefix = list(
    abbr="-p",
    help="Prefix for output files."
  )
)

scriptOptionalArgs <- list(
  score_column = list(
    default="logFC",
    help="Title to be put over the first subfigure."
  ),
  plot_title = list(
    default="Top gene sets",
    help="Title to be put over the first subfigure."
  ),
  msig_species = list(
    default="Mus musculus",
    help="Organism name where to look for gene symbols."
  ),
  msig_category = list(
    default="C5",
    help="Main division of MSigDB."
  ),
  msig_subcategory = list(
    default="BP",
    help="Subdivision of MSigDB."
  ),
  pAdjustMethod = list(
    default="none",
    help="Change this to BH for Bonferroni correction."
  ),
  pvalueCutoff = list(
    default=0.05,
    help="Cut-off for p-values."
  ),
  qvalueCutoff = list(
    default=1,
    help="An FDR of 0.05 is usually too stringent."
  ),
  minGSSize = list(
    default=10,
    help="Minimum number of genes in gene set."
  ),
  maxGSSize = list(
    default=500,
    help="Maximum number of genes in gene set."
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

for (pk in c("tidyr", "dplyr", "purr", "ggplot2", "cowplot", "msigdbr", "clusterProfiler", "rlang")){
  if(!(pk %in% (.packages()))){
    library(pk, character.only=TRUE)
  }
}

### Define a main function that will only be executed if called from command line
main <- function(opt){
  
  #' The main function of the script, executed only if called from command line.
  #' Calls subfunctions according to supplied command line arguments.
  #' 
  #' @param opt list. a named list of all command line options; will be passed on 
  #' 
  #' @return Not intended to return anything, but rather save outputs to files.
  
  outFile <- paste0(opt$outPrefix, opt$outFile) %>%
    gsub("/", "___", .)
  opt$outFile <- NULL
  opt$outPrefix <- NULL
  opt$commandRpath <- NULL
  opt$help <- NULL

  opt$geneSet <- download_ontologies(opt$msig_species, opt$msig_category, opt$msig_subcategory)
  if(opt$verbose){
    cat(paste0("Downloaded ", opt$msig_category, "/", opt$msig_subcategory, " for ", opt$msig_species, "\n"))
  }
  p <- do.call(plot_gsea, opt[!(names(opt) %in% c("msig_category", "msig_subcategory", "msig_species"))])
  
  if(opt$verbose){
    cat("Saving figure\n")
  }
  fig2pdf(p, outFile, height=9.6, width=7.2)
}

plot_gsea <- function(
  scoreTables,
  geneSet=NULL,
  emptyRes=NULL,
  score_column=NULL,
  universe=NULL,
  verbose=TRUE,
  ...
  ){
  
  #' Create a one-page figure showing top enriched gene sets (pathways) based on GSEA.
  #' 
  #' @description Downoads gene set information from MSigDB for a given species, runs
  #' Gene Set Enrichment Analysis and compiles a figure with 2 subplots: dotplot of gene
  #' sets and ridgeplot showing distribution of gene expression change in top gene sets.
  #' 
  #' @param scoreTables list. A named list of dataframes with gene ID in the first column,
  #' symbol in the second
  #' @param geneSet dataframe. Gene set membership of genes.
  #' @param emptyRes result object. An empty result object to be extended with enriched sets.
  #' @param score_column string. Name of column with scores. Third column if not specified.
  #' @param universe character vector. All genes in the organism.
  #' @param verbose logical. Whether progress messages should be printed.
  #' @param ... ellipse. Arguments to be passed on to the enricher function.
  #' @usage plot_gsea(scoreTables, geneSet=NULL, emptyRes=NULL, verbose=TRUE, ...)
  #' @return ggplot
  #' @details Dotplot of gene sets shows top enriched gene sets. Color corresponds to
  #' significance, while size shows...
  
  #' @examples
  #' plot_gsea(scoreTables)
  #' plot_gsea(scoreTables, geneSet)
  #' plot_gsea(scoreTables, geneSet, emptyRes)
  #' plot_gsea(scoreTables, geneSet, score_column="logFC", verbose=TRUE, pAdjustMethod="BH")
  
  if(is.null(geneSet)){
    if(verbose){
      cat("Downloading deafault ontology set\n")
    }
    geneSet <- download_ontologies()
  }

  if(is.null(universe)){
    universe <- unique(geneSet$gene_symbol)
  }
  
  if(is.null(emptyRes)){
    if(verbose){
      cat("Creating empty compareCluster object\n")
    }
    emptyRes <- create_empty_result_object()
  }
  
  if(verbose){
    cat("Looking for gene set enrichments\n")
  }
  enrichments <- scoreTables %>%
    map2(names(scoreTables), gsea_enrichment, geneSet, score_column=NULL, ...)

  if(verbose){
    cat("Plotting dotplot of top gene sets\n")
  }
  p1 <- enrichments %>%
    gsea_enrichments(geneSet, emptyRes, score_column, ...) %>%
    gsea_enrichdot()
  
  if(verbose){
    cat("Plotting ridgeplot of top gene sets\n")
  }
  p2 <- gsea_ridges(enrichments)

  if(verbose){
    cat("Combining subplots and saving figure\n")
  }
  p <- plot_grid(p1, p2, nrow=2, labels="AUTO")
  
  return(p)
}


gsea_enrichments <- function(enrichmentList, geneSet, emptyRes, score_column=NULL, ...){
  
  #' Do GSEA analysis for multiple conditions with ClusterProfiler and compare these
  #' results on a common dotplot.
  #' 
  #' @description Takes a list of tables. Names will appear as condition labels.
  #' ...
  #' 
  #' @param enrichmentList named list. Enrichment raults.
  #' @param geneSet dataframe. Gene set membership of genes.
  #' @param emptyRes result object. An empty result object to be extended with enriched sets.
  #' @param score_column string. Name of column with scores. Third column if not specified.
  #' @param ... ellipse. Arguments to be passed on to ClusterProfiler::GSEA.
  #' @usage single_gsea_enrichment(hitGenes, geneSet, emptyRes, score_column, ...)
  #' @return enrichment result

  #' @examples
  #' gsea_enrichments(scoreTable, geneSet, emptyRes, ...)
  #' gsea_enrichments(hitGenes, geneSet, emptyRes, score_column="logFC", pAdjustMethod="BH")


  compRes <- duplicate(emptyRes)
  compRes@compareClusterResult <- enrichmentList %>%
    map2(names(enrichmentList), gsea_enrichment, geneSet, score_column=NULL, ...) %>%
    map(function(x){x@result}) %>%
    bind_rows() %>%
    mutate(
      #Count = as.numeric(unlist(strsplit(GeneRatio, '/'))) / as.numeric(unlist(strsplit(BgRatio, '/'))))[seq(1, length(GeneRatio)*2, 2)],
      #Count = Count * 100,
      Description = gsub('GO_', '', Description),
      Description = gsub('_', ' ', Description)
    ) %>%
    arrange(p.adjust)

  return(compRes)
}


gsea_enrichment <- function(scoreTable, conditionName, geneSet, score_column=NULL, ...){
  
  #' Do GSEA analysis for a set of genes with numerical scores (e.g.: expression) with
  #' ClusterProfiler.
  #' 
  #' @description Takes a table with scores associated to gene names or symbols. 
  #' This score is typically logFC, but can also be p-value. GSEA is carried out after
  #' sorting based on score. 
  #' 
  #' @param scoreTable dataframe. A table with ID in the first column, symbol in the second.
  #' @param conditionName string. Name of the condition to be shown on plot.
  #' @param geneSet dataframe. Gene set membership of genes.
  #' @param score_column string. Name of column with scores. Third column if not specified.
  #' @param ... ellipse. Arguments to be passed on to ClusterProfiler::GSEA.
  #' @usage single_gsea_enrichment(hitGenes, geneSet, ...)
  #' @return encrichement result

  #' @examples
  #' gsea_enrichment(scoreTable, conditionName, geneSet)
  #' gsea_enrichment(hitGenes, conditionName, geneSet, score_column="logFC", pAdjustMethod="BH")

  cn <- colnames(scoreTable)
  if (is.null(score_column)){
    score_column <- cn[3]
  }
  browser()
  genesym <- cn[2]
  scoreTable <- scoreTable %>%
    distinct(across(one_of(c(genesym))), .keep_all = TRUE)
  hitGenes <- scoreTable[[score_column]]
  names(hitGenes) <- scoreTable[[hitGenes]]
  hitGenes <- hitGenes[order(hitGenes, decreasing=TRUE)]

  clusterProfiler::GSEA(hitGenes, TERM2GENE=geneSet, ...)
  enrichment@result <- enrichment@result %>%
    head(n=20) %>%
    mutate(
      Cluster = conditionName,
      group = conditionName
    )
    
  return(enrichment)
}


gsea_enrichdot <- function(enrichment, plot_title=opt$plot_title){

  #' Create a dotplot showing top enriched gene sets (pathways) for multiple hitlists.
  #' 
  #' @description A ClusterProfiler enrichment result for multiple hitlists is visualized
  #' side-by-side, as dotplot. 
  #' 
  #' @param enrichment. Result of clusterProfiler::enricher.
  #' @param plot_title string. Title of the figure.
  #' @usage single_hitlist_enrichdot(enrichment, plot_title="Top gene sets")
  #' @return ggplot
  #' @details Color corresponds to significance, while size shows gene count in hit list.
  
  #' @examples
  #' single_hitlist_enrichdot(enrichment)
  #' single_hitlist_enrichdot(enrichment, plot_title="Top gene sets")

resulTab <- enrichment@compareClusterResult
topsets <- rownames(resulTab[!duplicated(resulTab$ID), ])[1:30]
hitnames <- levels(resulTab$group)
clusterProfiler::dotplot(enrichment, showCategory=30) +
  labs(title=plot_title) +
  scale_color_gradientn(
    colors=rev(c('#2b8cbe', 'grey', '#e38071', '#e34a33', '#e31e00')),
    breaks=c(0.05, 0.01, 0.001, 0.0001),
    limits=c(0.00001, 1), trans='log10', oob = scales::squish
  ) +
  theme(
    axis.text.x=element_text(angle=30, hjust=1),
    axis.text.y=element_text(size=8)
  )
}

gsea_ridge <- function(enrichment, conditionName){
  
  #' Create a ridgeplot showing distribution of gene expression changes in top gene sets.
  #' 
  #' @description The top gene sets are scanned for hit genes, until 25 hit genes are
  #' collected. Membership of these top genes is shown in top gene sets.
  #' 
  #' @param enrichment ClusterProfiler result object. Result of an enrichment analysis.
  #' @param conditionName string. Name of the condition to be shown on plot.
  #' @usage gsea_ridge(enrichment, conditionName)
  #' @return ggplot
  #' @details ....
 
  #' @examples
  #' gsea_ridge(enrichment, conditionName)

  clusterProfiler::ridgeplot(enrichment, fill='pvalue') +
    xlab(conditionName) +
    scale_color_gradientn(
      colors=rev(c('#2b8cbe', 'grey', '#e38071', '#e34a33', '#e31e00')),
      breaks=c(0.05, 0.01, 0.001, 0.0001),
      limits=c(0.00001, 1), trans='log10', oob = scales::squish
    ) +
    scale_y_discrete(name="", limits=rev(topsets), labels=rev(topsets)) +
    scale_x_discrete(name="", labels=hitnames) +
    theme(
      axis.text.x=element_text(angle=30, hjust=1),
      axis.text.y=element_text(size=8)
    )
}

gsea_ridges <- function(enrichments){
  
  #' Create ridgeplots showing distribution of gene expression changes in top gene sets
  #' in all separate conditions.
  #' 
  #' @description ...
  #' 
  #' @param enrichments list. ClusterProfiler result objects.
  #' @return ggplot
  #' @details ....
 
  #' @examples
  #' gsea_ridge(enrichment, conditionName)

  grid_kws = map2(enrichments, names(enrichments), gsea_ridge)
  grid_kws[[nrow]] <- 2
  do.call(plot_grid, grid_kws))
  
}

download_ontologies <- function(msig_species=opt$msig_species, msig_category=opt$msig_category, msig_subcategory=opt$msig_subcategory){
  
  #' Download gene ontologies: a knowledge set is downloaded from MSigDB.
  #' 
  #' @description Downloeds gene ontologies (primarily the "Biological Process" section)
  #' as a dataframe, prettifies ontology names and removes extra columns not needed for
  #' the enrichment analysis.
  #' 
  #' @param msig_species string. Species for mapping gene symbols. 
  #' @param msig_category string. The section of MSigDB. 
  #' @param msig_subcategory string. A subsection of MSigDB. 
  #' @usage download_ontologies(msig_species="Mus musculus", msig_category="C5", msig_subcategory="BP")
  #' @return data frame with gene ontology membership of gene symbols

  #' @examples
  #' download_ontologies(msig_species="Mus musculus", msig_category="C5", msig_subcategory="BP")
  #' download_ontologies()

  geneSet <- msigdbr(species=msig_species, category=msig_category, subcategory=msig_subcategory)
  geneSet$gs_name <- gsub('GO_', '', geneSet$gs_name)
  geneSet$gs_name <- gsub('_', ' ', geneSet$gs_name)
  geneSet <- geneSet[,c('gs_name', 'gene_symbol')]
  return(geneSet)
}

create_empty_result_object <- function(){
  
  #' Create a ClusterProfile result object that will be extended with multi-enrichement
  #' result of our actual comparison.
  #' 
  #' @description This "empty" mock object is needed because the dotplot method only
  #' accepts dedicated objects.
  #' 
  #' @usage create_empty_result_object()
  #' @return ClusterProfile result object

  #' @examples
  #' create_empty_result_object()
  
  mydf <- data.frame(Entrez=c('1', '100', '1000', '100101467','100127206', '100128071'), group = c('A', 'A', 'A', 'B', 'B', 'B'))
  emptyRes <- compareCluster(Entrez~group, data=mydf, fun="enrichGO", 'org.Hs.eg.db')
  return(emptyRes)
}

# Ensuring command line connectivity by sourcing an argument parser
source(opt$commandRpath, local=TRUE)