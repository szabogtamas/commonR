# A script that takes orered or unordered hitlists and shows
# top gene ontologies, biological function or pathways enriched.

library(optparse)
library(docstring)
library(msigdbr)
library(clusterProfiler)
library(rlang)


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


# Check if mandatory arguments are present
if ( is.null(opt$inputfile) ) { 
  if ( opt$verbose ) { 
    write("Sorry, cannot proceed without a data table. Please provide a path to hitlists.\n", stderr())
  }
  checkpass <- FALSE
} else {
  checkpass <- TRUE
}

# Define helper functions

plotTopEnrichments <- function(
  compLists
)
{
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
  
  geneSet <-msigdbr(species="Homo sapiens", category="C5", subcategory="BP")
    compNames <- c('Array Cancer Up', 'Array Cancer Down')
    compLists <- list(rownames(tops[tops$adj.P.Val < 0.05 & tops$logFC > 0 ,]), rownames(tops[tops$adj.P.Val < 0.05 & tops$logFC < 0 ,]))
    richRes <- data.frame(Cluster=c(), group=c(), ID=c(), Description=c(), GeneRatio=c(), BgRatio=c(), pvalue=c(), p.adjust=c(), qvalue=c(), geneID=c(), Count=c())
    for (i in 1:length(compNames)){
    compname <- compNames[[i]]
    hitGenes <- unique(compLists[[i]])
    compname <- paste0(compname, ' (', length(hitGenes), ')')
    compNames[[i]] <- compname
    enrichment <- enricher(hitGenes, TERM2GENE=geneSet[,c('gs_name', 'gene_symbol')], universe=bgrdGenes)
    erik <- head(enrichment@result, n=25)
    erik$Cluster <- rep(compname, nrow(erik))
    erik$group <- rep(compname, nrow(erik))
    richRes <- rbind(erik, richRes)
    }
    richRes <- richRes[order(richRes$p.adjust), ]
    richRes$Count <- (as.numeric(unlist(strsplit(richRes$GeneRatio, '/'))) / as.numeric(unlist(strsplit(richRes$BgRatio, '/'))))[seq(1, length(richRes$GeneRatio)*2, 2)]
    richRes$Count <- richRes$Count * 100
    richRes$Description <- gsub('GO_', '', richRes$Description)
    richRes$Description <- gsub('NEGATIVE_REGULATION_OF_PATHWAY_RESTRICTED_', '', richRes$Description)
    richRes$Description <- gsub('ESTABLISHMENT_OF_PROTEIN_', '', richRes$Description)
    richRes$Description <- gsub('_', ' ', richRes$Description)
    compRes <- duplicate(emptyRes)
    richRes$group <- factor(richRes$group, levels=compNames)
    richRes$Cluster <- factor(richRes$Cluster, levels=compNames)
    compRes@compareClusterResult <- richRes
    richPlot <- clusterProfiler::dotplot(compRes, showCategory=10, by='count') + labs(title='Gene ontologies linked to hit genes')
    richPlot <- richPlot + scale_color_gradientn(colors=rev(c('#2b8cbe', 'grey', '#e38071', '#e34a33', '#e31e00')), breaks=c(0.05, 0.01, 0.001, 0.0001), limits=c(0.00001, 1), trans='log10', oob = scales::squish) + theme(axis.text.x=element_text(angle=30, hjust=1)) + scale_size_area(name='Percent\nin gene set')

  return(richPlot)
}

# Execute main function if mandatory arguments are supplied via the command line (otherwise print help message)
if ( checkpass ) { 
  clinicals <- read.table(opt$inputfile, header=TRUE, sep="\t")
  pdf(file=opt$km_out,7.2, 5.4, onefile=FALSE)
  print(plotTopEnrichments())
  if ( opt$verbose ) { 
    cat(paste0(Plots saved, "\n")) 
  }
  dev.off()
} else {
  print_help(parser)
}