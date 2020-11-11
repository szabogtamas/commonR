#!/usr/bin/env Rscript

scriptDescription <- "A script that takes a DE result with numerical changes under a condition and runs GSEA."

scriptMandatoryArgs <- list(
  scoreTables = list(
    abbr="-i",
    type="tables",
    help="Tables with unique ID and name as first two column and numeric as rest.",
    readoptions=list(stringsAsFactors=FALSE, sep="\t")
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
    default="BH",
    help="Change this to BH for Bonferroni correction."
  ),
  pvalueCutoff = list(
    default=1,
    help="Cut-off for p-values."
  ),
  minGSSize = list(
    default=10,
    help="Minimum number of genes in gene set."
  ),
  maxGSSize = list(
    default=500,
    help="Maximum number of genes in gene set."
  ),
  n_to_show = list(
    default=30,
    help="Number of top gene sets to show on dotplot."
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

for (pk in c("tidyr", "dplyr", "purrr", "ggplot2", "cowplot", "msigdbr", "clusterProfiler", "rlang")){
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
  #TODO: Sort out correct usage of ellipses
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

  additional_args <- list(...)
  plot_title <- additional_args[["plot_title"]]
  n_to_show <- additional_args[["n_to_show"]]
  additional_args <- additional_args[!(names(additional_args) %in% c("plot_title", "n_to_show"))]
  additional_args[["geneSet"]] <- geneSet
  additional_args[["score_column"]] <- score_column

  enrichments <- list()
  for (condition in names(scoreTables)){
    additional_args[["scoreTable"]] <- scoreTables[[condition]]
    additional_args[["conditionName"]] <- condition
    enrichments[[condition]] <- do.call(gsea_enrichments, additional_args)
  }
  
  if(verbose){
    cat("Plotting dotplot of top gene sets\n")
  }
  p1 <-  enrichments %>%
    merge_gsea_enrichments(emptyRes) %>%
    gsea_enrichdot(plot_title, n_to_show)

  if(verbose){
    cat("Plotting ridgeplot of top gene sets\n")
  }
  p2 <- gsea_ridges(enrichments, n_to_show)

  if(verbose){
    cat("Combining subplots and saving figure\n")
  }
  p <- plot_grid(p1, p2, nrow=2, labels="AUTO")
  
  invisible(p)

}


gsea_enrichments <- function(scoreTable, conditionName, geneSet, score_column=NULL, ...){
  
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
  genesym <- cn[2]
  scoreTable <- scoreTable %>%
    distinct(across(one_of(c(genesym))), .keep_all = TRUE)
  hitGenes <- scoreTable[[score_column]]
  names(hitGenes) <- scoreTable[[genesym]]
  hitGenes <- hitGenes[order(hitGenes, decreasing=TRUE)]
  
  enrichment <- clusterProfiler::GSEA(hitGenes, TERM2GENE=geneSet, ...)
  enrichment@result <- enrichment@result %>%
    mutate(
      absNES = abs(NES),
      Cluster = conditionName,
      group = conditionName
    ) %>%
    arrange(desc(absNES)) %>%
    head(n=20)
    
  return(enrichment)
}


merge_gsea_enrichments <- function(enrichmentList, emptyRes){
  
  #' Combine multiple GSEA enrichment results into a common object to be passed to
  #' dotplot visualization.
  #' 
  #' @description Takes a list of GSEA enrichment results, adds columns required by 
  #' dotplot and combines the result tables into a single table inside a result object.
  #' 
  #' @param enrichmentList named list. Enrichment results.
  #' @param emptyRes result object. An empty result object to be extended with enriched sets.
  #' @usage merge_gsea_enrichments(enrichmentList, emptyRes)
  #' @return enrichment result as a multi-condition result object

  #' @examples
  #' merge_gsea_enrichments(enrichmentList, emptyRes)


  compRes <- duplicate(emptyRes)
  enrichments <- enrichmentList %>%
    map(function(x) x@result) %>%
    bind_rows() %>%
    arrange(desc(absNES)) %>%
    mutate(
      Count = -log10(p.adjust),
      GeneRatio = -log10(p.adjust),
      BgRatio = "1/1"
    )
  compRes@compareClusterResult <- enrichments
  return(compRes)
}


gsea_enrichdot <- function(enrichment, plot_title="", n_to_show=30){

  #' Create a dotplot showing top enriched gene sets (pathways) for multiple hitlists.
  #' 
  #' @description A ClusterProfiler enrichment result for multiple conditions and after
  #' GSEA is visualized side-by-side, as dotplot. 
  #' 
  #' @param enrichment. Result of clusterProfiler::GSEA for multiple conditions.
  #' @param plot_title string. Title of the figure.
  #' @param n_to_show. Maximum number of enriched gene sets to show.
  #' @usage gsea_enrichdot(enrichment, plot_title="Top gene sets", n_to_show=30)
  #' @return ggplot
  #' @details Color corresponds to Normalized Enrichment score, while size shows
  #' significance.
  
  #' @examples
  #' gsea_enrichdot(enrichment)
  #' gsea_enrichdot(enrichment, plot_title="Top gene sets", n_to_show=30)

  clusterProfiler::dotplot(enrichment, showCategory=n_to_show, color="NES") +
    labs(
      title=plot_title,
      size="-log(p)"
    ) +
    scale_color_gradientn(
      colors=c('#2b8cbe', '#00bfff', 'grey', '#e38071', '#e31e00'),
      breaks=c(-2, 1, 0, 1, 2),
      limits=c(-2, 2), oob = scales::squish
    ) +
    theme(
      axis.text.x=element_text(angle=30, hjust=1),
      axis.text.y=element_text(size=8)
    )
}

gsea_ridge <- function(enrichment, conditionName, topsets){
  
  #' Create a ridgeplot showing distribution of gene expression changes in top gene sets.
  #' 
  #' @description ....
  #' 
  #' @param enrichment ClusterProfiler result object. Result of an enrichment analysis.
  #' @param conditionName string. Name of the condition to be shown on plot.
  #' @usage gsea_ridge(enrichment, conditionName)
  #' @return ggplot
  #' @details ....
 
  #' @examples
  #' gsea_ridge(enrichment, conditionName)
  
  enrichment@geneSets %>%
    .[topsets$ID] %>%
    enframe() %>%
    unnest() %>%
    mutate(
      gex = enrichment@geneList[value],
      condition = conditionName
    )
}


gsea_ridge2 <- function(enrichment, conditionName){
  
  #' Create a ridgeplot showing distribution of gene expression changes in top gene sets.
  #' 
  #' @description ....
  #' 
  #' @param enrichment ClusterProfiler result object. Result of an enrichment analysis.
  #' @param conditionName string. Name of the condition to be shown on plot.
  #' @usage gsea_ridge(enrichment, conditionName)
  #' @return ggplot
  #' @details ....
 
  #' @examples
  #' gsea_ridge(enrichment, conditionName)

  print(conditionName)
  save(enrichment, file="sdef.Rdata")
  clusterProfiler::ridgeplot(enrichment, fill="NES") +
    xlab(conditionName) +
    scale_color_gradientn(
      colors=rev(c("#2b8cbe", "grey", "#e38071", "#e34a33", "#e31e00")),
      breaks=c(0.05, 0.01, 0.001, 0.0001),
      limits=c(0.00001, 1), trans="log10", oob = scales::squish
    ) +
    scale_y_discrete(name="", limits=rev(topsets), labels=rev(topsets)) +
    scale_x_discrete(name="", labels=hitnames) +
    theme(
      axis.text.x=element_text(angle=30, hjust=1),
      axis.text.y=element_text(size=8)
    )
}

gsea_ridges <- function(enrichments, n_to_show=30){
  
  #' Create ridgeplots showing distribution of gene expression changes in top gene sets
  #' in all separate conditions.
  #' 
  #' @description ...
  #' 
  #' @param enrichment. Result of clusterProfiler::GSEA for multiple conditions.
  #' @return ggplot
  #' @details ....
 
  #' @examples
  #' gsea_ridges(enrichment)

  print("Ridging")
  print(names(enrichments))

  topsets <- list()
  for (enrn in names(enrichments)){
    enrichment <- enrichments[[enrn]]
    topset <- enrichment@result %>%
      mutate(
        absNES = abs(NES),
        group = enrn
      ) %>%
      arrange(desc(absNES)) %>%
      head(n = 20)
    topsets[[enrn]] <- topset
  }
  topsets <- topsets %>%
    bind_rows()
  
  topsets <- topsets %>%
    bind_rows() %>%
    pivot_wider(values_from=absNES, names_from=group, values_fill=0) %>%
    rowwise(ID) %>% 
    mutate(meanNES = mean(c_across(where(is.numeric)))) %>%
    arrange(desc(meanNES)) %>%
    slice(1:n_to_show) %>%
    .$ID %>%
    unique()

  print(topsets)

  enrichments %>%
    map2(names(enrichments), gsea_ridge, topsets) %>%
    bind_rows() %>%
    ggplot(aes(x=gex, fill=condition)) + 
    geom_density()  +
    facet_grid(rows = "name", scales = "free", switch="both") +
    scale_x_continuous(limits=c(-2,2)) +
    theme(
      strip.text.y.left = element_text(angle = 0),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()
    ) + 
    labs(x="logFC", y="")

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