#!/usr/bin/env Rscript

# A script that takes a hitlist and shows top gene ontologies

scriptMandatoryArgs <- c("hitGenes", "outFile")

scriptOptionalArgs <- list(
  plot_title = list(
    default='Top gene sets',
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
  qvalueCutoff = list(
    default=1,
    help="An FDR of 0.05 is usually too stringent."
  )

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("docstring"))

library(tidyr)
library(ggplot2)
library(cowplot)
library(msigdbr)
library(clusterProfiler)




if (!interactive()) {
  # Parse command line options if not sourced
  parser <- OptionParser()
  parser <- add_option(parser, c("-v", "--verbose"), action="store_true", default=FALSE,
                      help="Print some progress messages to stdout.")
  parser <- add_option(parser, c("-q", "--quietly"), action="store_false", 
                      dest="verbose", help="Create figures quietly, without printing to stdout.")
  parser <- add_option(parser, c("-i", "--hitGenes"), default=NULL, 
                      help="Comma separated list of hit genes.",
                      metavar="hit_genes")
  parser <- add_option(parser, c("-o", "--outFile"), default=NULL, 
                      help="File name prefix where to save output plots.",
                      metavar="msig_category")
  parser <- add_option(parser, c("-o", "--outFile"), default=NULL, 
                      help="File name prefix where to save output plots.",
                      metavar="msig_subcategory")
  parser <- add_option(parser, c("-o", "--plot_title"), default=NULL, 
                      help="File name prefix where to save output plots.",
                      metavar="out_prefix")
  parser <- add_option(parser, c("-o", "--msig_species"), default=NULL, 
                      help="File name prefix where to save output plots.",
                      metavar="out_prefix")
  parser <- add_option(parser, c("-o", "--outFile"), default=NULL, 
                      help="File name prefix where to save output plots.",
                      metavar="out_prefix")
  parser <- add_option(parser, c("-o", "--outFile"), default=NULL, 
                      help="File name prefix where to save output plots.",
                      metavar="out_prefix")
  opt <- parse_args(parser)

  print(opt)
}

### Define a main function that will only be executed if called from command line
main <- function(){
  geneSet <-msigdbr(species=msig_species, category=msig_category, subcategory=msig_subcategory)
  geneSet$gs_name <- gsub('GO_', '', geneSet$gs_name)
  geneSet$gs_name <- gsub('_', ' ', geneSet$gs_name)
  geneSet <- geneSet[,c('gs_name', 'gene_symbol')]

  enrichment <- single_enrichment(hitGenes, geneSet, pAdjustMethod=pAdjustMethod, qvalueCutoff=qvalueCutoff)
  p1 <- single_enrichdot(enrichment, plot_title=plot_title)
  p2 <- single_genedot(enrichment)

  p <- plot_grid(p1, p2, nrow=2, labels="AUTO")

  pdf(outFile, height=9.6, width=7.2)
  print(p)
  dev.off()
}

single_enrichment <- function(hitGenes, geneSet, pAdjustMethod=pAdjustMethod, qvalueCutoff=qvalueCutoff){
  
  #' Create a dotplot showing top enriched genes sets (pathways).
  #' 
  #' @description ...
  #' 
  #' @param compLists dataframe. Data to plot 
  #' @usage ...
  #' @return ...
  #' @details ...
  #' @examples
  #' ...

  enricher(hitGenes, TERM2GENE=geneSet, pAdjustMethod=pAdjustMethod, qvalueCutoff=qvalueCutoff)
}

single_enrichdot <- function(enrichment, plot_title=plot_title){

  #' Create a dotplot showing top enriched genes sets (pathways).
  #' 
  #' @description ...
  #' 
  #' @param compLists dataframe. Data to plot 
  #' @usage ...
  #' @return ...
  #' @details ...
  #' @examples
  #' ...

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

single_genedot <- function(enrichment){
  
  #' Create a dotplot showing top enriched genes sets (pathways).
  #' 
  #' @description ...
  #' 
  #' @param compLists dataframe. Data to plot 
  #' @usage ...
  #' @return ...
  #' @details ...
  #' @examples
  #' ...

  topenr <- enrichment@result[1:30, c("ID", "p.adjust", "BgRatio", "geneID")]
    rownames(topenr) <- topenr$ID

    detailedsets <- c()
    detailedgenes <- c()
    for (geneset in rownames(topenr)){
      if(length(detailedgenes) < 25){
        detailedgenes <- unique(c(detailedgenes, unlist(strsplit(topenr[geneset, "geneID"], "/"))))
        detailedsets <- c(detailedsets, geneset)
      }
    }

  geneFuns <- enrichment@result %>%
    .[detailedsets, c("ID", "p.adjust", "BgRatio", "geneID")] %>%
    transform(score = seq(length(detailedsets), 1, -1)) %>%
    transform(BgRatio = sapply(BgRatio, function(x){unlist(strsplit(x, "/"))[[1]]})) %>%
    transform(GeneSetSize = as.numeric(BgRatio)) %>%
    transform(geneID = as.character(geneID)) %>%
    transform(geneID = strsplit(geneID, "/")) %>%
    unnest(geneID) %>%
    group_by(geneID) %>%
    add_tally(wt=score) %>%
    mutate(genePriority = n) %>%
    arrange(desc(n))

  ggplot(data=geneFuns, aes(geneID, ID)) +
    geom_point(aes(size=GeneSetSize, color=genePriority)) +
    scale_y_discrete(name="", limits=rev(detailedsets), labels=rev(sapply(detailedsets, substr, 1, 35))) +
    scale_x_discrete(name="", limits=detailedgenes, labels=detailedgenes) +
    theme(axis.text.x=element_text(size=8, angle=30, hjust=1), legend.position="bottom")
}

if (!interactive()) {
  
  # Check if mandatory arguments are present
  if ( is.null(opt$inputfile) ) { 
    if ( opt$verbose ) { 
      write("Sorry, cannot proceed without a data table. Please provide a path to hitlists.\n", stderr())
    }
    checkpass <- FALSE
  } else {
    checkpass <- TRUE
  }

  # Execute main function if mandatory arguments set (otherwise print help message)
  if ( checkpass ) { 
    main()
  } else {
    print_help(parser)
  }
  
}