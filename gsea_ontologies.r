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
    default="none",
    help="Change this to BH for Bonferroni correction."
  ),
  pvalueCutoff = list(
    default=0.9,
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
  geneset_dist_plot = list(
    default="box",
    help="If boxplot or ridgeplot should show distribution of changes in gene sets."
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

for (
  pk in c(
    "tidyr", "dplyr", "purrr", "tibble", "ggplot2", "cowplot", "msigdbr",
    "clusterProfiler", "rlang"
  )
){
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
  #' @param score_column string. Name of column with scores. Third column if not specified.
  #' @param universe character vector. All genes in the organism.
  #' @param verbose logical. Whether progress messages should be printed.
  #' @param ... ellipse. Arguments to be passed on to the enricher function.
  #' @usage plot_gsea(scoreTables, geneSet=NULL, verbose=TRUE, ...)
  #' @return ggplot
  #' @details Dotplot of gene sets shows top enriched gene sets. Color corresponds to
  #' significance, while size shows...
  
  #' @examples
  #' plot_gsea(scoreTables)
  #' plot_gsea(scoreTables, geneSet)
  #' plot_gsea(scoreTables, geneSet)
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
  
  if(verbose){
    cat("Looking for gene set enrichments\n")
  }

  additional_args <- list(...)
  plot_title <- additional_args[["plot_title"]]
  n_to_show <- additional_args[["n_to_show"]]
  geneset_dist_plot <- additional_args[["geneset_dist_plot"]]
  additional_args <- additional_args[!(names(additional_args) %in% c("plot_title", "n_to_show", "geneset_dist_plot"))]
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
  p1 <-  gsea_enrichdot(enrichments, plot_title, n_to_show)

  if(verbose){
    cat("Plotting ridgeplot of top gene sets\n")
  }
  if (geneset_dist_plot == "ridge") {
    p2 <- gsea_ridges(enrichments, n_to_show)
  } else {
    p2 <- gsea_boxes(enrichments, n_to_show)
  }

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
  hitGenes <- hitGenes[!is.na(hitGenes)]

  enrichment <- clusterProfiler::GSEA(hitGenes, TERM2GENE=geneSet, ...)
  enrichment@result <- enrichment@result %>%
    mutate(
      absNES = abs(NES),
      group = conditionName
    ) %>%
    arrange(desc(absNES)) %>%
    head(n=20)
  
  return(enrichment)
}


gsea_enrichdot <- function(enrichmentList, plot_title="", n_to_show=20){

  #' Create a dotplot showing top enriched gene sets (pathways) for multiple hitlists.
  #' 
  #' @description A ClusterProfiler enrichment result for multiple conditions and after
  #' GSEA is visualized side-by-side, as dotplot. 
  #' 
  #' @param enrichmentList named list. Enrichment results.
  #' @param plot_title string. Title of the figure.
  #' @param n_to_show. Maximum number of enriched gene sets to show.
  #' @usage gsea_enrichdot(enrichmentList, plot_title="Top gene sets", n_to_show=20)
  #' @return ggplot
  #' @details Color corresponds to Normalized Enrichment score, while size shows
  #' significance.
  
  #' @examples
  #' gsea_enrichdot(enrichmentList)
  #' gsea_enrichdot(enrichmentList, plot_title="Top gene sets")

  save(enrichmentList, file = "rich_list.RData")
  enrichmentList <- enrichmentList %>%
    map(function(x) x@result) %>%
    bind_rows()
    
  topsets <- enrichmentList %>%
    .$Description %>%
    unique() %>%
    head(n_to_show)

  rich_plot_dot <- enrichmentList %>%
    arrange(desc(absNES)) %>%
    filter(Description %in% topsets)

  save(rich_plot_dot, file = "rich_dot.RData")
  rich_plot_dot %>%
    ggplot(aes(x=group, y=Description, fill=NES, size=absNES)) +
    geom_dotplot(stackdir="center") +
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

gsea_ridge_rich <- function(enrichment, conditionName, topsets){
  
  #' Take top gene sets in a GSEA result and assign gene expression change (or other score)
  #' to each gebe set.
  #' 
  #' @description Collects scores for all genes in top gene sets for all conditions and
  #' merges this information into a dataframe that can later be used for plotting.
  #' 
  #' @param enrichment ClusterProfiler result object. Result of an enrichment analysis.
  #' @param conditionName string. Name of the condition to be shown on plot.
  #' @param topsets character vector. Vector specifying the most enriched gene sets.
  #' @usage gsea_ridge_rich(enrichment, conditionName, topsets)
  #' @return ggplot
 
  #' @examples
  #' gsea_ridge_rich(enrichment, conditionName, topsets)
  
  enrichment@geneSets %>%
    .[topsets] %>%
    enframe() %>%
    unnest() %>%
    mutate(
      gex = enrichment@geneList[value],
      condition = conditionName
    )
}

gsea_boxes <- function(enrichments, n_to_show=30){
  
  #' Create boxplots showing distribution of gene expression changes in top gene sets
  #' in all separate conditions.
  #' 
  #' @description Gene expression changes of all genes in a given gene set are shown on 
  #' a density plot. Color of histograms corresponds to condition.
  #' 
  #' @param enrichment. Result of clusterProfiler::GSEA for multiple conditions.
  #' @usage gsea_boxes(enrichments, n_to_show=30)
  #' @return ggplot
  
  #' @examples
  #' gsea_boxes(enrichment)

  topsets <- list()
  for (enrn in names(enrichments)){
    enrichment <- enrichments[[enrn]]
    topsets[[enrn]] <- enrichment@result %>%
      mutate(
        absNES = abs(NES),
        group = enrn
      ) %>%
      arrange(desc(absNES)) %>%
      head(n_to_show)
  }
  
  topsets <- topsets %>%
    bind_rows() %>%
    pivot_wider(values_from=absNES, names_from=group, values_fill=0) %>%
    rowwise(ID) %>% 
    mutate(meanNES = mean(c_across(where(is.numeric)))) %>%
    arrange(desc(meanNES)) %>%
    head(n_to_show) %>%
    .$ID %>%
    unique()

  rich_plot_box <- enrichments %>%
    map2(names(enrichments), gsea_ridge_rich, topsets) %>%
    bind_rows() %>%
    mutate(
      name = factor(name, levels=topsets),
      name = substr(name, 1, 35)
    ) 
  
  save(rich_plot_box, file = "rich_box.RData") # Save but tsv!
  rich_plot_box %>%
    ggplot(aes(x=name, y=gex, fill=condition)) + 
    geom_boxplot(position=position_dodge(1), outlier.shape = NA) +
    theme(
      axis.ticks = element_blank(),
      axis.text.x=element_text(size=7, angle=30, hjust=1),
      legend.position = "left"
    ) + 
    ylim(-1, 1) +
    labs(x="", y="logFC")

}

gsea_ridges <- function(enrichments, n_to_show=30){
  
  #' Create ridgeplots showing distribution of gene expression changes in top gene sets
  #' in all separate conditions.
  #' 
  #' @description Gene expression changes of all genes in a given gene set are shon on 
  #' a density plot. Color of histograms corresponds to condition.
  #' 
  #' @param enrichment. Result of clusterProfiler::GSEA for multiple conditions.
  #' @usage gsea_ridges(enrichments, n_to_show=30)
  #' @return ggplot
  
  #' @examples
  #' gsea_ridges(enrichment)

  topsets <- list()
  for (enrn in names(enrichments)){
    enrichment <- enrichments[[enrn]]
    topsets[[enrn]] <- enrichment@result %>%
      mutate(
        absNES = abs(NES),
        group = enrn
      ) %>%
      arrange(desc(absNES)) %>%
      head(n_to_show)
  }
  
  topsets <- topsets %>%
    bind_rows() %>%
    pivot_wider(values_from=absNES, names_from=group, values_fill=0) %>%
    rowwise(ID) %>% 
    mutate(meanNES = mean(c_across(where(is.numeric)))) %>%
    arrange(desc(meanNES)) %>%
    head(n_to_show) %>%
    .$ID %>%
    unique()

  enrichments %>%
    map2(names(enrichments), gsea_ridge_rich, topsets) %>%
    bind_rows() %>%
    mutate(
      name = factor(name, levels=topsets),
      name = substr(name, 1, 35)
    ) %>%
    ggplot(aes(x=gex, fill=condition)) + 
    geom_density()  +
    facet_grid(rows = "name", scales = "free", switch="both") +
    scale_x_continuous(limits=c(-2,2)) +
    theme(
      strip.text.y.left = element_text(size=8, angle = 0),
      strip.background = element_rect(colour="white", fill="#ffffff"),
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

# Ensuring command line connectivity by sourcing an argument parser
source(opt$commandRpath, local=TRUE)