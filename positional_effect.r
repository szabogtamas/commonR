#!/usr/bin/env Rscript

scriptDescription <- "A script that takes a DE result and check overrepresented genomic locations via GSEA."

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
  )
)

scriptOptionalArgs <- list(
  outPrefix = list(
    abbr="-p",
    default="",
    help="Prefix for output files."
  ),
  conditionOrder = list(
    default=NULL,
    type="vector",
    help="Order of experimental conditions on plots."
  ),
  conditionColors = list(
    default=NULL,
    type="vector",
    help="Colors for points belonging to a given condition on the FC distribution plot."
  ),
  score_column = list(
    default="logFC",
    help="Column that has the values to be shown."
  ),
  plot_title = list(
    default="Top genomic postions",
    help="Title to be put over the first subfigure."
  ),
  msig_species = list(
    default="Mus musculus",
    help="Organism name where to look for gene symbols."
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


#' The main function of the script, executed only if called from command line.
#' Calls subfunctions according to supplied command line arguments.
#' 
#' @param opt list. a named list of all command line options; will be passed on 
#' 
#' @return Not intended to return anything, but rather save outputs to files.
main <- function(opt){
  
  outFile <- paste0(opt$outPrefix, opt$outFile) %>%
    gsub("/", "___", .)
  opt$outFile <- NULL
  opt$outPrefix <- NULL
  opt$commandRpath <- NULL
  opt$help <- NULL

  opt$geneSet <- download_ontologies(opt$msig_species, "C1", NULL)
  p <- do.call(plot_gsea, opt[!(names(opt) %in% c("msig_species"))])
  
  if(opt$verbose){
    cat("Combining subplots and saving figure\n")
  }
  p <- plot_grid(p$genedot, p$setchange, nrow=2, labels="AUTO")

  fig2pdf(p, outFile, height=9.6, width=7.2)
}


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
#' @return list two output plots and also the enrichments
#' @details Dotplot of gene sets shows top enriched gene sets. Color corresponds to
#' significance, while size shows...

#' @examples
#' plot_gsea(scoreTables)
#' plot_gsea(scoreTables, geneSet)
#' plot_gsea(scoreTables, geneSet)
#' plot_gsea(scoreTables, geneSet, score_column="logFC", verbose=TRUE, pAdjustMethod="BH")
plot_gsea <- function(
  scoreTables,
  geneSet=NULL,
  score_column=NULL,
  universe=NULL,
  verbose=TRUE,
  ...
  ){

  if(is.null(geneSet)){
    if(verbose){
      cat("Downloading deafault ontology set\n")
    }
    geneSet <- download_ontologies()
  }

  if(is.null(universe)){
    universe <- unique(geneSet$gene_symbol)
  }

  if (!exists("conditionOrder")){
    conditionOrder <- names(scoreTables)
  }
  if(is.null(conditionOrder)){
    conditionOrder <- names(scoreTables)
  }

  if (!exists("conditionColors")){
    conditionColors <- default_colors
  }
  if(is.null(conditionColors)){
    conditionColors <- default_colors
  }
  
  if(verbose){
    cat("Looking for gene set enrichments\n")
  }

  additional_args <- list(...)
  plot_title <- additional_args[["plot_title"]]
  n_to_show <- additional_args[["n_to_show"]]
  conditionOrder <- additional_args[["conditionOrder"]]
  conditionColors <- additional_args[["conditionColors"]]
  geneset_dist_plot <- additional_args[["geneset_dist_plot"]]
  additional_args <- additional_args[!(names(additional_args) %in% c(
    "plot_title", "n_to_show", "geneset_dist_plot", "conditionOrder", "conditionColors"
  ))] 
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
  p1 <-  gsea_enrichdot(enrichments, plot_title, n_to_show, conditionOrder)
  
  if(verbose){
    cat("Plotting distibution of changes in top gene sets\n")
  }
  if (geneset_dist_plot == "ridge") {
    p2 <- gsea_ridges(enrichments, n_to_show, topsets=p1$topsets)
  } else {
    p2 <- gsea_boxes(enrichments, n_to_show, conditionOrder, conditionColors, topsets=p1$topsets)
  }
  
  invisible(list(genedot=p1$plot, setchange=p2, enrichments=enrichments, topsets=p1$topsets))

}


#' Show fold changes accross chromosomes
#' 
#' @description Takes a table with scores associated to gene names or symbols. 
#' This score is typically logFC, but can also be p-value. GSEA is carried out after
#' sorting based on score. 
#' 
#' @param scoreTable dataframe. A table with ID in the first column, symbol in the second.
#' @param conditionName string. Name of the condition to be shown on plot.
#' @param geneSet dataframe. Gene set membership of genes.
#' @param score_column string. Name of column with scores. Third column if not specified.
#' @usage gsea_local_enrichments(hitGenes, geneSet, score_column=NULL)
#' @return encrichement result

#' @examples
#' gsea_local_enrichments(scoreTable, conditionName, geneSet)
#' gsea_local_enrichments(hitGenes, conditionName, geneSet, score_column="logFC")
fc_local <- function(scoreTable, conditionName, geneSet, score_column=NULL){

gene_positions <- opt$geneSet %>%
  mutate(
    chromosome = stringr:::str_extract(gs_name, "chr(X|Y|\\d*)"),
    chromosomen = stringr:::str_replace(chromosome, "chrX", "23"),
    chromosomen = stringr:::str_replace(chromosomen, "chrY", "24"),
    chromosomen = stringr:::str_replace(chromosomen, "chr", ""),
    chromosomen = as.numeric(chromosomen),
    band = stringr:::str_replace(gs_name, "chr\\d*", ""),
    band = stringr:::str_replace_all(band, "q", "-"),
    band = stringr:::str_replace_all(band, "[a-zA-Z]", ""),
    band = as.numeric(band),
    band = ifelse(band < 0, band + 10, band - 10),
    location = 100 * chromosomen + band
  )

p <- scoreTables %>%
  imap(~mutate(.x, condition = .y)) %>%
  bind_rows() %>%
  left_join(gene_positions, by = c(Symbol = "gene_symbol")) %>%
  mutate(
    condition = factor(condition),
    chromosome = factor(chromosome)
  ) %>%
  ggplot(aes(x=band, y=logFC, color=logFC)) +
    geom_point() +
    facet_grid(rows = vars(condition), cols = vars(chromosome), scales = "free") +
    scale_color_gradientn(
      colors=c('#2b8cbe', '#00bfff', 'grey', 'grey', 'grey', '#e38071', '#e31e00'),
      breaks=c(-4, -2, -1, 0, 1, 2, 4),
      limits=c(-4, 4), oob = scales::squish
    )

  return(p)
}


#' Show fold changes accross chromosomes
#' 
#' @description Takes a table with scores associated to gene names or symbols. 
#' This score is typically logFC, but can also be p-value. GSEA is carried out after
#' sorting based on score. 
#' 
#' @param scoreTable dataframe. A table with ID in the first column, symbol in the second.
#' @param conditionName string. Name of the condition to be shown on plot.
#' @param geneSet dataframe. Gene set membership of genes.
#' @param score_column string. Name of column with scores. Third column if not specified.
#' @usage gsea_local_enrichments(hitGenes, geneSet, score_column=NULL)
#' @return encrichement result

#' @examples
#' gsea_local_enrichments(scoreTable, conditionName, geneSet)
#' gsea_local_enrichments(hitGenes, conditionName, geneSet, score_column="logFC")
plot_local_changes <- function(scoreTable, conditionName, geneSet, score_column=NULL){

gene_positions <- opt$geneSet %>%
  mutate(
    is_mitochondr = grepl("MT", gs_name),
    chromosome = str_extract(gs_name, "chr(X|Y|\\d*)"),
    chromosome = ifelse(is_mitochondr, "MT", chromosome),
    chromosomen = str_replace(chromosome, "chrX", "23"),
    chromosomen = str_replace(chromosomen, "chrY", "24"),
    chromosomen = str_replace(chromosomen, "chr", ""),
    chromosomen = as.numeric(chromosomen),
    band  = ifelse(is_mitochondr, "10", gs_name),
    band = str_replace(band, "chr\\d*", ""),
    band = str_replace_all(band, "q", "-"),
    band = str_replace_all(band, "[a-zA-Z]", ""),
    band = as.numeric(band),
    band = ifelse(band < 0, band + 10, band - 10),
    location = 100 * chromosomen + band
  )


scoreTables %>%
  imap(~mutate(.x, condition = .y)) %>%
  bind_rows() %>%
  left_join(gene_positions, by = c(Symbol = "gene_symbol")) %>%
  filter(!is.na(chromosome)) %>%
  mutate(
    condition = factor(condition),
    chromosome = factor(chromosome, levels = c(
      "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
      "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
      "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
      "chr22", "chrX", "chrY", "MT" 
    ))
  ) %>%
  group_by(condition, chromosome, band) %>%
  mutate(
    median_logFC = median(logFC),
    q1_logFC = quantile(logFC, 0.05),
    q3_logFC = quantile(logFC, 0.95)
  ) %>%
  ggplot(aes(x=band)) +
    facet_grid(rows = vars(condition), cols = vars(chromosome), scales = "free", switch="both") +
    geom_point(aes(y = logFC, color = logFC)) +
    geom_ribbon(aes(ymin = q1_logFC, ymax = q3_logFC), fill = "#F9E076") +
    geom_line(aes(y = median_logFC)) +
    scale_color_gradientn(
      colors=c('#2b8cbe', '#00bfff', 'grey', 'grey', 'grey', '#e38071', '#e31e00'),
      breaks=c(-4, -2, -1, 0, 1, 2, 4),
      limits=c(-4, 4), oob = scales::squish
    ) +
    theme_bw() +
    theme(
      strip.text.y.left = element_text(size=8, angle = 0),
      strip.text.x.bottom = element_text(size=6, angle = 60, margin=margin(t=2, r=7, b=7, l=7)),
      strip.background = element_rect(colour="white", fill="#ffffff"),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) + 
    labs(x="", y="logFC")

  return(p)
}


#' Do GSEA analysis for a set of genes with numerical scores (e.g.: expression) with
#' ClusterProfiler on chromosomal locations (bands).
#' 
#' @description Takes a table with scores associated to gene names or symbols. 
#' This score is typically logFC, but can also be p-value. GSEA is carried out after
#' sorting based on score. 
#' 
#' @param scoreTable dataframe. A table with ID in the first column, symbol in the second.
#' @param conditionName string. Name of the condition to be shown on plot.
#' @param geneSet dataframe. Gene set membership of genes.
#' @param score_column string. Name of column with scores. Third column if not specified.
#' @usage gsea_local_enrichments(hitGenes, geneSet, score_column=NULL)
#' @return encrichement result

#' @examples
#' gsea_local_enrichments(scoreTable, conditionName, geneSet)
#' gsea_local_enrichments(hitGenes, conditionName, geneSet, score_column="logFC")
gsea_local_enrichments <- function(scoreTable, conditionName, geneSet, score_column=NULL){

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

  enrichment <- clusterProfiler::GSEA(
    hitGenes,
    TERM2GENE=geneSet,
    pAdjustMethod="none",
    pvalueCutoff=1
  )
  enrichment@result <- enrichment@result %>%
    mutate(
      absNES = abs(NES),
      group = conditionName
    ) %>%
    arrange(desc(absNES))
  
  return(enrichment)
}


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
download_ontologies <- function(msig_species=opt$msig_species, msig_category=opt$msig_category, msig_subcategory=opt$msig_subcategory){

  geneSet <- msigdbr(species=msig_species, category=msig_category, subcategory=msig_subcategory)
  geneSet$gs_name <- gsub('GO_', '', geneSet$gs_name)
  geneSet$gs_name <- gsub('KEGG_', '', geneSet$gs_name)
  geneSet$gs_name <- gsub('HALLMARK_', '', geneSet$gs_name)
  geneSet$gs_name <- gsub('_', ' ', geneSet$gs_name)
  geneSet <- geneSet[,c('gs_name', 'gene_symbol')]
  return(geneSet)
}

# Ensuring command line connectivity by sourcing an argument parser
source(opt$commandRpath, local=TRUE)