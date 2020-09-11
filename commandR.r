#!/usr/bin/env Rscript

# Interface to command line based on optparse, but extending it with lists and dataframes

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(docstring))

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
      if ("readoptions" %in% names(rg) ) {
        rg[["readoptions"]] <- NULL
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
        if (rg[["type"]] %in% c("vector", "nested", "table") ) {
          if (rg[["type"]] == "vector") {
            opt[[rn]] <- unlist(strsplit(opt[[rn]], ",", fixed=TRUE))
          } else {
            if (rg[["type"]] == "nested") {
              nl <- list()
              for (x in unlist(strsplit(opt[[rn]], ":", fixed=TRUE))){
                x <- unlist(strsplit(x, ",", fixed=TRUE))
                nl[[x[1]]] <- x[2:length(x)]
              }
              opt[[rn]] <- nl
            } else {
              opt[[rn]] <- do.call(read.csv, c(list(opt[[rn]]), rg[["readoptions"]]))
            }
            
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
    main(opt)
  } else {
    print_help(parser)
  }
  
}