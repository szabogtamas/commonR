#!/usr/bin/env Rscript

scriptDescription <- "A script that takes a hitlist and shows top gene ontologies"

scriptMandatoryArgs <- list(
  hitGenes = list(
    abbr="-i",
    type="nested",
    help="A comma separated list of genes, usually a hitlist."
  ),
  outFile = list(
    abbr="-o",
    help="Prefix for output files."
  )
)

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

### Define a main function that will only be executed if called from command line
main <- function(verbose, ...){
  print(verbose)
  single_enrichment(...)
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

  print(hitGenes)
  print(paste0(hitGenes, "hi"))
}

if (!interactive()) {
  
  # Initialize parser with verbosity and description of script
  parser <- OptionParser(usage=paste0("%prog [options]\nDescription:\n  ", scriptDescription))
  parser <- add_option(
    parser,
    c("-v", "--verbose"),
    action="store_true",
    default=FALSE,
    help="Print some progress messages to stdout."
    )
  parser <- add_option(
    parser,
    c("-q", "--quietly"),
    action="store_false",
    dest="verbose",
    help="Create figures quietly, without printing to stdout."
    )

  # Add arguments to parser
  an <- 0
  for (al in list(scriptMandatoryArgs, scriptOptionalArgs)){
    for (rgn in names(al)){
      rg <- al[[rgn]]
      if (an < 1){
        rg[["default"]] <- NULL
      }

      rga <- paste0("--", rgn)
      if ("abbr" %in% names(rg) ) {
        rga <- c(rg[["abbr"]], rga)
        rg[["abbr"]] <- NULL
      }

      if ("type" %in% names(rg) ) {
        rg[["type"]] <- NULL
      }

      rl <- list(parser, rga)
      rl <- c(rl, rg)
      parser <- do.call(add_option, rl)
      an <- an +1
    }
  }

  # Parse command line options and split up lists or nested lists
  opt <- parse_args(parser)
  for (rn in names(c(scriptMandatoryArgs, scriptOptionalArgs))){
    rg <- c(scriptMandatoryArgs, scriptOptionalArgs)[[rn]]
    if ("type" %in% names(rg) ) {
        if (rg[["type"]] %in% c("vector", "nested") ) {
          if (rg[["type"]] == "vector") {
            opt[[rn]] <- unlist(strsplit(opt[[rn]], ",", fixed=TRUE))
          } else {
            opt[[rn]] <- unlist(strsplit(opt[[rn]], ":", fixed=TRUE))
          }
        }
      }
  }
  

  # Check if mandatory arguments are present
  passed_args <- opt[names(scriptMandatoryArgs)]
  if (any(is.na(names(passed_args)))) {
    if (opt$verbose) { 
      write("Sorry, cannot proceed without all mandatory arguments.\n", stderr())
    }
    checkpass <- FALSE
  } else {
    checkpass <- TRUE
  }

  # Execute main function if mandatory arguments are set (otherwise print help message)
  if (checkpass) { 
    do.call(main, opt)
  } else {
    print_help(parser)
  }
  
}