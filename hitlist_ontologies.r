#!/usr/bin/env Rscript

scriptDescription <- "A script that takes a hitlist and shows top gene ontologies"

scriptMandatoryArgs <- list(
  hitGenes = list(
    abbr="-i",
    type="list",
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

opt <- list()
for (rn in names(scriptOptionalArgs)){
  opt[[rn]] <- scriptOptionalArgs[[rn]][["default"]]
}


suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(docstring))

### Define a main function that will only be executed if called from command line
main <- function(verbose, outFile, msig_species, msig_category, msig_subcategory, ...){

  if(verbose){
    cat(paste0("Downloading ", msig_category, "/", msig_subcategory, "for ", msig_species))
  }
  geneSet <-msigdbr(species=msig_species, category=msig_category, subcategory=msig_subcategory)
  geneSet$gs_name <- gsub('GO_', '', geneSet$gs_name)
  geneSet$gs_name <- gsub('_', ' ', geneSet$gs_name)
  geneSet <- geneSet[,c('gs_name', 'gene_symbol')]

  if(verbose){
    cat("Looking for gene set enrichments")
  }
  enrichment <- single_enrichment(...)

  if(verbose){
    cat("Plotting dotplot of top gene sets")
  }
  p1 <- single_enrichdot(enrichment, ...)
  
  if(verbose){
    cat("Plotting dotplot of top genes")
  }
  p2 <- single_genedot(enrichment, ...)

  if(verbose){
    cat("Combining subplots and saving figure")
  }
  pdf(paste0(outFile, ".pdf"), height=9.6, width=7.2)
  plot_grid(p1, p2, nrow=2, labels="AUTO")
  dev.off()
}

single_enrichment <- function(hitGenes, geneSet, ...){
  
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
  enricher(hitGenes, TERM2GENE=geneSet, universe=universe, pAdjustMethod=pAdjustMethod, qvalueCutoff=qvalueCutoff, pvalueCutoff=pvalueCutoff, minGSSize=minGSSize, maxGSSize=maxGSSize)
}

single_enrichdot <- function(enrichment, plot_title=opt$plot_title){

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

plotMultipleEnrichments <- function(
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
            nl <- list()
            for (x in unlist(strsplit(opt[[rn]], ":", fixed=TRUE))){
              x <- unlist(strsplit(x, ",", fixed=TRUE))
              nl[[x[1]]] <- x[2:length(x)]
            }
            opt[[rn]] <- nl
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