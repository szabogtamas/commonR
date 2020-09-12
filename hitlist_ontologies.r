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
  pvalueCutoff = list(
    default=0.05,
    help="Cut-off for p-values."
  ),
  qvalueCutoff = list(
    default=1,
    help="An FDR of 0.05 is usually too stringent."
  ),
  minGSSize = list(
    default=10,
    help="Minimum number of genes in gene set."
  ),
  maxGSSize = list(
    default=500,
    help="Maximum number of genes in gene set."
  )
)

opt <- list()
for (rn in names(scriptOptionalArgs)){
  opt[[rn]] <- scriptOptionalArgs[[rn]][["default"]]
}

for (pk in c("tidyr", "dplyr", "ggplot2", "cowplot", "msigdbr", "clusterProfiler", "rlang")){
  if(!(pk %in% (.packages()))){
    library(pk, character.only=TRUE)
  }
}

### Define a main function that will only be executed if called from command line
main <- function(opt){
  
  outFile <- opt$outFile
  opt$outFile <- NULL
  opt$help <- NULL

  opt$geneSet <- download_ontologies(opt$msig_species, opt$msig_category, opt$msig_subcategory)
  if(opt$verbose){
    cat(paste0("Downloaded ", opt$msig_category, "/", opt$msig_subcategory, " for ", opt$msig_specie, "\n"))
  }

  if (length(opt$hitGenes) > 1){
    if(opt$verbose){
      cat("Multiple hitlists supplied. Plotting comparison on common plot\n")
    }
    p <- do.call(plot_enrichment_for_multiple_hitlist, opt[!(names(opt) %in% c("plot_title", "msig_category", "msig_subcategory", "msig_species"))])
  } else {
    if(opt$verbose){
      cat("Single hitlists supplied\n")
    }
    p <- do.call(plot_enrichment_for_single_hitlist, opt[!(names(opt) %in% c("plot_title", "msig_category", "msig_subcategory", "msig_species"))])
  }

  cat("Saving figure\n")
  pdf(paste0(outFile, ".pdf"), height=9.6, width=7.2)
  print(p)
  dev.off()
}

plot_enrichment_for_single_hitlist <- function(hitGenes, geneSet=NULL, verbose=TRUE, ...){

  if(is.null(geneSet)){
    if(verbose){
      cat("Downloading deafault ontology set\n")
    }
    geneSet <- download_ontologies()
  }

  if(verbose){
    cat("Looking for gene set enrichments\n")
  }
  hitGenes <- unlist(hitGenes)
  enrichment <- single_enrichment(hitGenes, geneSet, ...)

  if(verbose){
    cat("Plotting dotplot of top gene sets\n")
  }
  p1 <- single_enrichdot(enrichment)
  
  if(verbose){
    cat("Plotting dotplot of top genes\n")
  }
  p2 <- single_genedot(enrichment)

  if(verbose){
    cat("Combining subplots and saving figure\n")
  }
  p <- plot_grid(p1, p2, nrow=2, labels="AUTO")
  return(p)
}

plot_enrichment_for_multiple_hitlist <- function(hitGenes, geneSet=NULL, emptyRes=NULL, verbose=TRUE, ...){
  
  if(is.null(geneSet)){
    if(verbose){
      cat("Downloading deafault ontology set\n")
    }
    geneSet <- download_ontologies()
  }
  
  if(is.null(emptyRes)){
    if(verbose){
      cat("Creating empty compareCluster object\n")
    }
    emptyRes <- create_empty_result_object()
  }
  richRes <- data.frame(Cluster=c(), group=c(), ID=c(), Description=c(), GeneRatio=c(), BgRatio=c(), pvalue=c(), p.adjust=c(), qvalue=c(), geneID=c(), Count=c())
  
  if(verbose){
    cat("Looking for gene set enrichments\n")
  }
  compNames <- list()
  i <- 0
  for (compname in names(hitGenes)){
    i <- i + 1
    hitGenIt <- unique(hitGenes[[compname]])
    compname <- paste0(compname, ' (', length(hitGenIt), ')')
    compNames[[i]] <- compname
    enrichment <- single_enrichment(hitGenes, geneSet, ...)
    erik <- head(enrichment@result, n=20)
    erik$Cluster <- rep(compname, nrow(erik))
    erik$group <- rep(compname, nrow(erik))
    richRes <- rbind(erik, richRes)
    }
  richRes <- richRes[order(richRes$p.adjust), ]
  richRes$Count <- (as.numeric(unlist(strsplit(richRes$GeneRatio, '/'))) / as.numeric(unlist(strsplit(richRes$BgRatio, '/'))))[seq(1, length(richRes$GeneRatio)*2, 2)]
  richRes$Count <- richRes$Count * 100
  richRes$Description <- gsub('GO_', '', richRes$Description)
  richRes$Description <- gsub('_', ' ', richRes$Description)
  compRes <- duplicate(emptyRes)
  richRes$group <- factor(richRes$group, levels=compNames)
  richRes$Cluster <- factor(richRes$Cluster, levels=compNames)
  compRes@compareClusterResult <- richRes
  return(compRes)
}

download_ontologies <- function(msig_species=opt$msig_species, msig_category=opt$msig_category, msig_subcategory=opt$msig_subcategory){
  geneSet <- msigdbr(species=msig_species, category=msig_category, subcategory=msig_subcategory)
  geneSet$gs_name <- gsub('GO_', '', geneSet$gs_name)
  geneSet$gs_name <- gsub('_', ' ', geneSet$gs_name)
  geneSet <- geneSet[,c('gs_name', 'gene_symbol')]
  return(geneSet)
}

create_empty_result_object <- function(){
  mydf <- data.frame(Entrez=c('1', '100', '1000', '100101467','100127206', '100128071'), group = c('A', 'A', 'A', 'B', 'B', 'B'))
  emptyRes <- compareCluster(Entrez~group, data=mydf, fun="enrichGO", 'org.Hs.eg.db')
  return(emptyRes)
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

  universe <- unique(geneSet$gene_symbol)
  enricher(hitGenes, TERM2GENE=geneSet, universe=universe, ...)
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
  #' @description The top gene sets are scanned for hit genes, until 25 hit genes are collected. Membership of these top genes is shown in top gene sets.
  #' 
  #' @param enrichment ClusterProfiler result object. Result of an enrichment analysis.
  #' @usage single_genedot(enrichment)
  #' @return ggplot
  #' @details Dot size shows how many genes the gene set consists of. Color shows how gene sets linked to a certain gene ranked.
  #' @examples
  #' single_genedot(enrichment)

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

source("/home/szabo/dev_packages/commonR/commandR.r")