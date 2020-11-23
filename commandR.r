#!/usr/bin/env Rscript

# Interface to command line based on optparse, but extending it with lists and dataframes

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(docstring))

default_colors <- c(
  '#1a476f', '#90353b', '#55752f', '#e37e00', '#6e8e84', '#c10534',
  '#938dd2', '#cac27e', '#a0522d', '#7b92a8', '#2d6d66', '#9c8847',
  '#bfa19c', '#ffd200', '#d9e6eb'
)

### Patching Pheatmap with diagonal column labels;
### Idea from https://www.thetopsites.net/article/54919955.shtml
draw_colnames_30 <- function (coln, gaps, ...) {
    coord = pheatmap:::find_coordinates(length(coln), gaps)
    x = coord$coord - 0.5 * coord$size
    res = grid:::textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 30, gp = grid:::gpar(...))
    return(res)}

assignInNamespace(x="draw_colnames", value="draw_colnames_30", ns=asNamespace("pheatmap"))


#' Save a plot, preferably one page, into pdf.
#' 
#' @param figure plot (ggplot or cowplot). The plot to be saved.
#' @param filename_base string. File name, with path and prefix, but no extension. 
#' @param height integer. Height of the output canvas.
#' @param width integer. Width of the output canvas. 
#' 
#' @return NULL.
fig2pdf <- function(figure, filename_base, height, width){
  pdf(paste0(filename_base, ".pdf"), height=height, width=width)
  print(figure)
  dev.off()
  invisible(NULL)
}


#' Save a dataframe to tsv. Rownames bacome first columns. If rownames are not human-
#' friendly IDs, the second column can be a more readable mapping, provided by relabels.
#' 
#' @param tab dta.fram. The table to be saved.
#' @param filename_base string. File name, with path and prefix, but no extension. 
#' @param primary_id string. Column header for row names.
#' @param secondary_id string. Column header for second column with human-friendly labels.
#' @param relabels named vector. Mapping from IDs to human-friendly labels. 
#' 
#' @return NULL.
genetab2tsv <- function(tab, filename_base, primary_id="GeneID", secondary_id="Symbol", relabels=NULL){
  original_cols <- colnames(tab)
  if(!is.null(relabels)){
    tab[[secondary_id]] <- relabels[rownames(tab)]
    original_cols <- c(secondary_id, original_cols)
  }
  original_cols <- c(primary_id, original_cols)
  tab %>% 
    rownames_to_column(var=primary_id) %>% 
    select(original_cols) %>% 
    write.table(
      paste0(filename_base, ".tsv"),
      quote=FALSE,
      sep="\t",
      row.names=FALSE
    )
  invisible(NULL)
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
    if ("type" %in% names(rg) & !is.null(opt[[rn]])) {
        if (rg[["type"]] %in% c("vector", "nested", "table", "tables") ) {
          if (rg[["type"]] == "vector") {
            opt[[rn]] <- unlist(strsplit(opt[[rn]], ",", fixed=TRUE))
          } else {
            if (rg[["type"]] == "nested") {
              nl <- list()
              sl <- unlist(strsplit(opt[[rn]], ":", fixed=TRUE))
              for (x in sl){
                x <- unlist(strsplit(x, ",", fixed=TRUE))
                if (length(sl) > 1){
                  nl[[x[1]]] <- x[2:length(x)]
                } else {
                  nl[[1]] <- x
                }
              }
              opt[[rn]] <- nl
            } else {
              if (rg[["type"]] == "tables") {
                nl <- list()
                sl <- unlist(strsplit(opt[[rn]], ",", fixed=TRUE))
                n <- 0
                for (x in sl){
                  n <- n+1
                  x <- unlist(strsplit(x, ":", fixed=TRUE))
                  if (length(x) > 1){
                    nl[[x[1]]] <- do.call(read.csv, c(list(x[2]), rg[["readoptions"]]))
                  } else {
                    nl[[paste0("Condition_", n)]] <- do.call(read.csv, c(list(x[1]), rg[["readoptions"]]))
                  }
                }
                opt[[rn]] <- nl
              } else {
                opt[[rn]] <- do.call(read.csv, c(list(opt[[rn]]), rg[["readoptions"]]))
              }
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