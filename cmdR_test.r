#!/usr/bin/env Rscript

scriptDescription <- "A script that takes a hitlist and shows top gene ontologies"

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
)

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(docstring))

if (!interactive()) {
  # Parse command line options if not sourced
  parser <- OptionParser(usage=paste0("Purpose: ", scriptDescription, "\nUsage: %prog [options]")
  parser <- add_option(parser, c("-v", "--verbose"), action="store_true", default=FALSE,
                      help="Print some progress messages to stdout.")
  parser <- add_option(parser, c("-q", "--quietly"), action="store_false", 
                      dest="verbose", help="Create figures quietly, without printing to stdout.")
  parser <- add_option(parser, c("-i", "--hitGenes"), default=NULL, 
                      help="Comma separated list of hit genes.",
                      metavar="hit_genes")
  for (rgn in names(scriptOptionalArgs)){
    print(rgn)
    rg <- scriptOptionalArgs[[rgn]]
    rl <- list(
      parser,
      paste0("--", deparse(substitute(rgn)))
    )
    rl <- c(rl, rg)
    parser <- do.call(add_option, rl)
  }

  print(parser)
}

### Define a main function that will only be executed if called from command line
main <- function(hitGenes, ...){
  single_enrichment(hitGenes)
}

single_enrichment <- function(hitGenes, ...){
  
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

  print(paste0(hitGenes, "hi"))
}

if (!interactive()) {
  
  # Check if mandatory arguments are present
  print("hi")
  if ( is.null(opt$hitGenes) ) { 
    if ( opt$verbose ) { 
      write("Sorry, cannot proceed without a data table. Please provide a path to hitlists.\n", stderr())
    }
    checkpass <- FALSE
  } else {
    checkpass <- TRUE
  }

  print(checkpass)
  # Execute main function if mandatory arguments set (otherwise print help message)
  if ( checkpass ) { 
    do.call(main, opt)
  } else {
    print_help(parser)
  }
  
}