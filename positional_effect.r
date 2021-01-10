#!/usr/bin/env Rscript

not_called_by_another <- FALSE
source("gsea_ontologies.r", local=TRUE)

scriptDescription <- "A script that takes a DE result and checks overrepresented genomic locations via GSEA."

scriptMandatoryArgs <- list(
  scoreTables = list(
    abbr="-i",
    type="tables",
    help="Tables with unique ID and name as first two columns and numeric as rest.",
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
    "tidyr", "dplyr", "stringr", "purrr", "tibble", "ggplot2", "cowplot", "msigdbr",
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
  p <- do.call(plot_positional, opt)
  
  if(opt$verbose){
    cat("Combining subplots and saving figure\n")
  }
  p <- plot_grid(p$genedot, p$topology, nrow=2, labels="AUTO")

  fig2pdf(p, outFile, height=9.6, width=7.2)
}


#' Create a one-page figure showing topological distribution of gene expression and the
#' probability of certain chromosomal regions being enriched based on GSEA.
#' 
#' @description Downoads gene set information from MSigDB for a given species, runs
#' Gene Set Enrichment Analysis and compiles a figure with 2 subplots: dotplot of gene
#' sets and a Manhattan-like plot showing topological distribution of gene expression.
#' 
#' @param scoreTables list. A named list of dataframes with gene ID in the first column,
#' symbol in the second and a score; typically logFC.
#' @param geneSet dataframe. Gene set membership of genes.
#' @param msig_species string. Species of origin for genes.
#' @param score_column string. Name of column with scores. The logFC column if not specified.
#' @param universe character vector. All genes in the organism.
#' @param verbose logical. Whether progress messages should be printed.
#' @param pAdjustMethod string. If adjustment of GSEA p-value is needed, the method can be specfied.
#' @param pvalueCutoff float. The cutoff value for enrichments.
#' @param n_to_show. Maximum number of enriched gene sets to show.
#' @param ... ellipse. Arguments to be passed on to the enricher function.
#' @usage plot_positional(scoreTables, geneSet=NULL, verbose=TRUE, pAdjustMethod="none", pvalueCutoff=0.05, ...)
#' @return list. Two output plots and also the enrichments.
#' @details Dotplot of gene sets shows top enriched gene sets. Color corresponds to
#' significance, while size shows the direction of the change.

#' @examples
#' plot_positional(scoreTables)
#' plot_positional(scoreTables, geneSet)
#' plot_positional(scoreTables, geneSet)
#' plot_positional(scoreTables, geneSet, score_column="logFC", verbose=TRUE, pAdjustMethod="BH")
plot_positional <- function(
  scoreTables,
  geneSet=NULL,
  msig_species="Mus musculus",
  score_column=NULL,
  universe=NULL,
  verbose=TRUE,
  pAdjustMethod="none",
  pvalueCutoff=0.05,
  n_to_show=20,
  ...
  ){

  if(is.null(geneSet)){
    if(verbose){
      cat("Downloading deafault ontology set\n")
    }
    geneSet <- download_ontologies(msig_species, "C1", NULL)
  }

  if(is.null(universe)){
    universe <- unique(geneSet$gene_symbol)
  }

  additional_args <- list(...)
  plot_title <- additional_args[["plot_title"]]
  conditionOrder <- additional_args[["conditionOrder"]]
  conditionColors <- additional_args[["conditionColors"]]
 
  if(is.null(conditionOrder)){
    conditionOrder <- names(scoreTables)
  }
  if(is.null(conditionColors)){
    conditionColors <- default_colors
  }
  
  if(verbose){
    cat("Enrichments of chromosomal locations\n")
  }
  additional_args <- additional_args[!(names(additional_args) %in% c(
    "plot_title", "conditionOrder", "conditionColors"
  ))] 
  additional_args[["geneSet"]] <- geneSet
  additional_args[["score_column"]] <- score_column
  additional_args[["pAdjustMethod"]] <- pAdjustMethod
  additional_args[["pvalueCutoff"]] <- pvalueCutoff

  enrichments <- list()
  for (condition in names(scoreTables)){
    additional_args[["scoreTable"]] <- scoreTables[[condition]]
    additional_args[["conditionName"]] <- condition
    enrichments[[condition]] <- do.call(gsea_enrichments, additional_args)
  }
  
  if(verbose){
    cat("Plotting dotplot of top gene sets\n")
  }
  p1 <-  gsea_enrichdot(enrichments, plot_title, n_to_show, conditionOrder, gs_hjust=0)
  
  if(verbose){
    cat("Plotting topological distibution of fold changes\n")
  }
  p2 <- plot_local_changes(scoreTables, geneSet, score_column=score_column)
  
  invisible(list(genedot=p1$plot, topology=p2, enrichments=enrichments, topsets=p1$topsets))

}


#' Show fold changes accross chromosomes
#' 
#' @description Takes a list of tables with scores of genes for different conditions and
#' Shows gene expression changes along chromosomal positions to help spot structural effects. 
#' This score is typically logFC.
#' 
#' @param scoreTables dataframe. A list of tables with scores of genes for different conditions.
#' @param geneSet dataframe. Gene set membership of genes.
#' @param score_column string. Name of column with scores. The logFC column if not specified.
#' @usage plot_local_changes(scoreTables, geneSet, score_column=NULL)
#' @return ggplot showing fold change of genes on a band

#' @examples
#' plot_local_changes(scoreTables, geneSet)
#' plot_local_changes(scoreTables, geneSet, score_column="logFC")
plot_local_changes <- function(scoreTables, geneSet, score_column=NULL){
  
  if(is.null(score_column)){
    score_column <- "logFC"
  }
  if(score_column == "logFC"){
    colrenames <- setNames(c(score_column), c("logFC"))
  } else {
    colrenames <- setNames(c(score_column, "logFC"), c("logFC", "logFC1"))
  }
  
gene_positions <- geneSet %>%
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
  rename(
    any_of(colrenames)
  ) %>%
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
    geom_point(aes(y = logFC, color = logFC), size=0.5) +
    geom_ribbon(aes(ymin = q1_logFC, ymax = q3_logFC), fill = "#F9E076") +
    geom_line(aes(y = median_logFC)) +
    scale_color_gradientn(
      colors=c('#2b8cbe', '#00bfff', 'grey', 'grey', 'grey', '#e38071', '#e31e00'),
      breaks=c(-4, -2, -1, 0, 1, 2, 4),
      limits=c(-4, 4), oob = scales::squish
    ) +
    theme_bw() +
    theme(
      strip.text.y.left = element_text(size=8, angle=0, vjust=0),
      strip.text.x.bottom = element_text(size=6, angle=60, margin=margin(t=2, r=7, b=7, l=7)),
      strip.background = element_rect(colour="white", fill="#ffffff"),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.1, "lines")
    ) + 
    labs(x="", y="", color=score_column)
}

# Ensuring command line connectivity by sourcing an argument parser
not_called_by_another <- TRUE
source(opt$commandRpath, local=TRUE)