#!/usr/bin/env Rscript

# A script that takes a hitlist and shows top gene ontologies

library(optparse)
library(docstring)
library(tidyr)
library(ggplot2)
library(cowplot)
library(msigdbr)
library(clusterProfiler)

# Parse command line options
parser <- OptionParser()
parser <- add_option(parser, c("-v", "--verbose"), action="store_true", default=FALSE,
                     help="Print some progress messages to stdout.")
parser <- add_option(parser, c("-q", "--quietly"), action="store_false", 
                     dest="verbose", help="Create figures quietly, without printing to stdout.")
parser <- add_option(parser, c("-f", "--inputfile"), default=NULL, 
                     help="Data file containing gene names, optionally scores and labels.",
                     metavar="path_to_hitlist")
opt <- parse_args(parser)

#TODO: Add more options, design input format (list? dataframe?)


# Check if mandatory arguments are present
if ( is.null(opt$inputfile) ) { 
  if ( opt$verbose ) { 
    write("Sorry, cannot proceed without a data table. Please provide a path to hitlists.\n", stderr())
  }
  checkpass <- FALSE
} else {
  checkpass <- TRUE
}


args = commandArgs(trailingOnly=TRUE)

outFile <- args[1]
hitGenes <- unlist(strsplit(args[2], ",", fixed=TRUE))

plot_title <- 'Top gene sets'
msig_species <- "Mus musculus"
msig_category <- "C5"
msig_subcategory <- "BP"

pAdjustMethod <- "none"
qvalueCutoff <- 1

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
  enricher(hitGenes, TERM2GENE=geneSet, pAdjustMethod=pAdjustMethod, qvalueCutoff=qvalueCutoff)
}

single_enrichdot <- function(enrichment, plot_title=plot_title){
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
  main()
}


# Execute main function if mandatory arguments are supplied via the command line (otherwise print help message)
if ( checkpass ) { 
  clinicals <- read.table(opt$inputfile, header=TRUE, sep="\t")
  pdf(file=opt$km_out,7.2, 5.4, onefile=FALSE)
  print(plotTopEnrichments()) # TODO: add option for returning a singe pdf file, separate figures or the figures as objects ("pickled?")
  if ( opt$verbose ) { 
    cat(paste0(Plots saved, "\n")) 
  }
  dev.off()
} else {
  print_help(parser)
}