#!/usr/bin/env Rscript

scriptDescription <- "A script that shows DE fold changes on a KEGG pathway map."

scriptMandatoryArgs <- list(
  scoreTable = list(
    abbr="-i",
    type="table",
    help="Table with ENSEMBL IDs and fold change values.",
    readoptions=list(stringsAsFactors=FALSE, sep="\t")
  ),
  keggPath = list(
    abbr="-p",
    help="KEGG ID of pathwas (e.g. hsa04110).",
    readoptions=list(stringsAsFactors=FALSE, sep="\t")
  ),
  outFile = list(
    abbr="-o",
    help="Base name for output files."
  )
)

scriptOptionalArgs <- list(
  score_column = list(
    default="logFC",
    help="Column that has the values to be shown."
  ),
  plot_title = list(
    default="Gene expression changes",
    help="Title to be put over the first subfigure."
  ),
  species = list(
    default="hsa",
    help="Species to be used when mapping gene Ids."
  ),
  scorecol = list(
    default="logFC",
    help="Column name containing scores (fold changes)."
  ),
  gene_idcol = list(
    default="GeneID",
    help="Column name containing gene IDs."
  ),
  gene_mapping = list(
    default="ENSEMBL",
    help="Gene ID type."
  ),
  commandRpath = list(
    default="commandR.r",
    help="Path to command line connectivity script (if not in cwd)."
  ),
  out_path = list(
    default=NULL,
    help="Path to command line connectivity script (if not in cwd)."
  )
)

opt <- list()
for (rn in names(scriptOptionalArgs)){
  opt[[rn]] <- scriptOptionalArgs[[rn]][["default"]]
}

for (
  pk in c(
    "tidyr", "dplyr", "pathview", "ggplot2", "cowplot"
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

  opt$outFile <- gsub("/", "___", opt$outFile)
  outFile <- opt$outFile
  plot_title <- opt$plot_title
  opt$plot_title <- NULL
  opt$commandRpath <- NULL
  opt$help <- NULL

  if(!is.null(opt$out_path)){
    setwd(opt$out_path)
    opt$out_path <- NULL
  }

  if(opt$verbose){
    cat(paste0("Drawing pathway map for ", keggPath, " (", outFile, ")\n"))
  }
  do.call(draw_kegg_path, opt)

  if(opt$verbose){
    cat(paste0("Adding title to plot\n"))
  }
  p <- ggplot() +
    theme_void() +
    draw_image(paste0(keggPath, ".", outFile, ".png"))
        
  fig2pdf(p, outFile, height=9.6, width=7.2)
}


#' Draw gene expression changes on a KEGG pathway map.
#' 
#' @description Draw gene expression changes on a KEGG pathway map using pathview.
#' 
#' @param scoreTable data.frame. A dataframe with two mandatory columns: gene IDs and scores.
#' @param keggPath string. KEGG ID of the pathway to be shown.
#' @param species string. KEGG species string (e.g. hsa).
#' @param gene_idcol string. Name of column that features gene IDs.
#' @param scorecol string. Name of column that features expression changes.
#' @param gene_mapping string. What kind of gene IDs are given.
#' @param outFile string. suffix to be added in output file name.
#' @param ... ellipse. Arguments to be passed on to pathview.
#' @usage draw_kegg_path(scoreTable, keggPath, species, gene_idcol, scorecol, gene_mapping, outFile, ...)
#' @return Does not return the plot, just generates an image in the working directory.
#' @details Uses the neat pathway representation provided by KEGG to zoom in to DE genes
#' in a given pathway. A major limitation is that it generates a png file (pdf is not as
#' neat) and customization of the plot is hard.

#' @examples
#' draw_kegg_path(scoreTable, keggPath, species, gene_idcol, scorecol, gene_mapping, outFile)
#' draw_kegg_path(scoreTable, keggPath, species, gene_idcol, scorecol, gene_mapping, outFile, kegg.native = FALSE)
draw_kegg_path <- function(
  scoreTable,
  keggPath,
  species,
  gene_idcol,
  scorecol,
  gene_mapping,
  outFile,
  ...
  ){

    if(gene_mapping == "ENSEMBL"){
        scoreTable$GeneID <- scoreTable[[gene_idcol]]
        scoreTable <- scoreTable %>%
            mutate(
                GeneID = lapply(GeneID, function(x) strsplit(x, ".", fixed = TRUE)),
                GeneID = lapply(GeneID, function(x) unlist(x)[[1]])
            ) %>%
            unnest()
        gene_idcol 
    }

    custargs <- list(...)
    custargs[["gene.data"]] <- setNames(scoreTable[[scorecol]], scoreTable[[gene_idcol]])
    custargs[["pathway.id"]] <- keggPath
    custargs[["species"]] <- species
    custargs[["gene.idtype"]] <- gene_mapping
    custargs[["out.suffix"]] <- outFile

    do.call(pathview, custargs)

}


# Ensuring command line connectivity by sourcing an argument parser
source(opt$commandRpath, local=TRUE)