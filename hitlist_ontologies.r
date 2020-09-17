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

for (pk in c("tidyr", "dplyr", "ggplot2", "ggplotify", "cowplot", "msigdbr", "clusterProfiler", "rlang")){
  if(!(pk %in% (.packages()))){
    library(pk, character.only=TRUE)
  }
}

### Define a main function that will only be executed if called from command line
main <- function(opt){
  
  #' The main function of the script, executed only if called from command line.
  #' Calls subfunctions according to supplied command line arguments.
  #' 
  #' @param opt list. a named list of all command line options; will be passed on 
  #' 
  #' @return Not intended to return enything, but rather save outputs to files.
  
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

plot_enrichment_for_single_hitlist <- function(hitGenes, geneSet=NULL, universe=NULL, verbose=TRUE, ...){
  
  #' Create a one-page figure showing top enriched genes sets (pathways) for a single gene set.
  #' 
  #' @description Downoads gene set information from MSigDB for a given species, runs
  #' overrepresentation analysis and compiles a figure with 2 subplots: dotplot of gene
  #' sets and dotplot showing gene set membership for top (best known) genes.
  #' 
  #' @param hitGenes character vector. A vector of gene symbols to be queried.
  #' @param geneSet dataframe. Gene set membership of genes.
  #' @param verbose logical. Whether progress messages should be printed.
  #' @param ... ellipse. Arguments to be passed on to the enricher function.
  #' @usage plot_enrichment_for_single_hitlist(hitGenes, geneSet=NULL, verbose=TRUE, ...)
  #' @return ggplot
  #' @details Dotplot of gene sets show top enriched gene sets. Color corresponds to
  #' significance, while size shows gene count in hit list. Dotplot of gene set mebership
  #' features dot where genes belong to a given gene set. Size corresponds to gene set size.
  
  #' @examples
  #' plot_enrichment_for_single_hitlist(hitGenes)
  #' plot_enrichment_for_single_hitlist(hitGenes, geneSet)
  #' plot_enrichment_for_single_hitlist(hitGenes, geneSet, verbose=FALSE)
  #' plot_enrichment_for_single_hitlist(hitGenes, geneSet, verbose=FALSE, pAdjustMethod="BH")

  if(is.null(geneSet)){
    if(verbose){
      cat("Downloading deafault ontology set\n")
    }
    geneSet <- download_ontologies()
  }
  
  if(is.null(universe)){
    universe <- unique(geneSet$gene_symbol)
  }

  if(verbose){
    cat("Looking for gene set enrichments\n")
  }
  hitGenes <- unlist(hitGenes)
  enrichment <- single_hitlist_enrichment(hitGenes, geneSet, universe, ...)

  if(verbose){
    cat("Plotting dotplot of top gene sets\n")
  }
  p1 <- single_hitlist_enrichdot(enrichment)
  
  if(verbose){
    cat("Plotting dotplot of top genes\n")
  }
  p2 <- single_hitlist_genedot(enrichment)

  if(verbose){
    cat("Combining subplots and saving figure\n")
  }
  p <- plot_grid(p1, p2, nrow=2, labels="AUTO")
  return(p)
}

plot_enrichment_for_multiple_hitlist <- function(hitGenes, geneSet=NULL, emptyRes=NULL, universe=NULL, verbose=TRUE, ...){
  
  #' Create a one-page figure showing top enriched gene sets (pathways) for mulitple gene sets.
  #' 
  #' @description Downoads gene set information from MSigDB for a given species, runs
  #' overrepresentation analysis and compiles a figure with 2 subplots: dotplot of gene
  #' sets and dotplot showing gene set membership for top (best known) genes.
  #' 
  #' @param hitGenes character vector. A vector of gene symbols to be queried.
  #' @param geneSet dataframe. Gene set membership of genes.
  #' @param emptyRes result object. An empty result object to be extended with enriched sets.
  #' @param universe character vector. All genes in the organism.
  #' @param verbose logical. Whether progress messages should be printed.
  #' @param ... ellipse. Arguments to be passed on to the enricher function.
  #' @usage plot_enrichment_for_single_hitlist(hitGenes, geneSet=NULL, emptyRes=NULL, verbose=TRUE, ...)
  #' @return ggplot
  #' @details Dotplot of gene sets show top enriched gene sets. Color corresponds to
  #' significance, while size shows gene count in hit list. Dotplot of gene set mebership
  #' features dot where genes belong to a given gene set. Size corresponds to gene set size.
  
  #' @examples
  #' plot_enrichment_for_single_hitlist(hitGenes)
  #' plot_enrichment_for_single_hitlist(hitGenes, geneSet)
  #' plot_enrichment_for_single_hitlist(hitGenes, geneSet, emptyRes)
  #' plot_enrichment_for_single_hitlist(hitGenes, geneSet, verbose=TRUE, pAdjustMethod="BH")
  
  if(is.null(geneSet)){
    if(verbose){
      cat("Downloading deafault ontology set\n")
    }
    geneSet <- download_ontologies()
  }

  if(is.null(universe)){
    universe <- unique(geneSet$gene_symbol)
  }
  
  if(is.null(emptyRes)){
    if(verbose){
      cat("Creating empty compareCluster object\n")
    }
    emptyRes <- create_empty_result_object()
  }
  
  if(verbose){
    cat("Looking for gene set enrichments\n")
  }
  enrichment <- multi_hitlist_enrichment(hitGenes, geneSet, emptyRes, universe, ...)

  if(verbose){
    cat("Plotting dotplot of top gene sets\n")
  }
  p1 <- multi_hitlist_enrichdot(enrichment)
  
  if(verbose){
    cat("Plotting dotplot of top genes\n")
  }
  p2 <- multi_hitlist_genedot(enrichment)

  if(verbose){
    cat("Combining subplots and saving figure\n")
  }
  p <- plot_grid(p1, p2, nrow=2, labels="AUTO")
  return(p)
}

download_ontologies <- function(msig_species=opt$msig_species, msig_category=opt$msig_category, msig_subcategory=opt$msig_subcategory){
  
  #' Download gene ontologies: a knowledge set is downloaded from MSigDB.
  #' 
  #' @description Downloeds gene ontologies (primarily the "Biological Process" section)
  #' as a dataframe, prettifies ontology names and removes extra columns not needed for
  #' the enrichment analysis.
  #' 
  #' @param msig_species string. Species for mapping gene symbols. 
  #' @param msig_category string. The section of MSigDB. 
  #' @param msig_subcategory string. A subsection of MSigDB. 
  #' @usage download_ontologies(msig_species="Mus musculus", msig_category="C5", msig_subcategory="BP")
  #' @return data frame with gene ontology membership of gene symbols

  #' @examples
  #' download_ontologies(msig_species="Mus musculus", msig_category="C5", msig_subcategory="BP")
  #' download_ontologies()

  geneSet <- msigdbr(species=msig_species, category=msig_category, subcategory=msig_subcategory)
  geneSet$gs_name <- gsub('GO_', '', geneSet$gs_name)
  geneSet$gs_name <- gsub('_', ' ', geneSet$gs_name)
  geneSet <- geneSet[,c('gs_name', 'gene_symbol')]
  return(geneSet)
}

create_empty_result_object <- function(){
  
  #' Create a ClusterProfile result object that will be extended with multi-enrichement
  #' result of our actual comparison.
  #' 
  #' @description This "empty" mock object is needed because the dotplot method only
  #' accepts dedicated objects.
  #' 
  #' @usage create_empty_result_object()
  #' @return ClusterProfile result object

  #' @examples
  #' create_empty_result_object()
  
  mydf <- data.frame(Entrez=c('1', '100', '1000', '100101467','100127206', '100128071'), group = c('A', 'A', 'A', 'B', 'B', 'B'))
  emptyRes <- compareCluster(Entrez~group, data=mydf, fun="enrichGO", 'org.Hs.eg.db')
  return(emptyRes)
}

single_hitlist_enrichment <- function(hitGenes, geneSet, ...){
  
  #' Do overrepresentation analysis for a single gene list with ClusterProfiler.
  #' 
  #' @description Needs a gene list and a knowledge set of gene ontology memberships to
  #' to do ORA. Extra arguments will be passed on to clusterProfiler::enricher
  #' 
  #' @param hitGenes character vector. A vector of gene symbols to be queried.
  #' @param geneSet dataframe. Gene set membership of genes.
  #' @param ... ellipse. Arguments to be passed on to ClusterProfiler::enricher.
  #' @usage single_enrichment(hitGenes, geneSet, ...)
  #' @return enrichment result

  #' @examples
  #' single_enrichment(hitGenes, geneSet, ...)
  #' single_enrichment(hitGenes, geneSet, pAdjustMethod="BH")

  clusterProfiler::enricher(hitGenes, TERM2GENE=geneSet, ...)
}

multi_hitlist_enrichment <- function(hitGenes, geneSet, emptyRes, ...){
  
  #' Do overrepresentation analysis for multiple gene lists with ClusterProfiler.
  #' 
  #' @description Takes a nested gene list and a knowledge set of gene ontology
  #' memberships to to do ORA.
  #' 
  #' @param hitGenes character vector. A vector of gene symbols to be queried.
  #' @param geneSet dataframe. Gene set membership of genes.
  #' @param emptyRes result object. An empty result object to be extended with enriched sets.
  #' @param ... ellipse. Arguments to be passed on to ClusterProfiler::enricher.
  #' @usage single_enrichment(hitGenes, geneSet, ...)
  #' @return enrichment result

  #' @examples
  #' single_enrichment(hitGenes, geneSet, emptyRes)
  #' single_enrichment(hitGenes, geneSet, emptyRes, pAdjustMethod="BH")

  richRes <- data.frame(
    Cluster=c(),
    group=c(),
    ID=c(),
    Description=c(),
    GeneRatio=c(),
    BgRatio=c(),
    pvalue=c(),
    p.adjust=c(),
    qvalue=c(),
    geneID=c(),
    Count=c()
    )

  compNames <- list()
  i <- 0
  for (compname in names(hitGenes)){
    i <- i + 1
    hitGenIt <- unique(hitGenes[[compname]])
    compname <- paste0(compname, ' (', length(hitGenIt), ')')
    compNames[[i]] <- compname
    enrichment <- single_hitlist_enrichment(hitGenIt, geneSet, ...)
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

single_hitlist_enrichdot <- function(enrichment, plot_title=opt$plot_title){

  #' Create a dotplot showing top enriched gene sets (pathways).
  #' 
  #' @description A ClusterProfiler enrichment result is visualized as dotplot. 
  #' 
  #' @param enrichment. Result of  clusterProfiler::enricher.
  #' @param plot_title string. Title of the figure.
  #' @usage single_hitlist_enrichdot(enrichment, plot_title="Top gene sets")
  #' @return ggplot
  #' @details Color corresponds to significance, while size shows gene count in hit list.
  
  #' @examples
  #' single_hitlist_enrichdot(enrichment)
  #' single_hitlist_enrichdot(enrichment, plot_title="Top gene sets")

topsets <- rownames(enrichment@result)[1:30]
clusterProfiler::dotplot(enrichment, showCategory=30) +
  labs(title=plot_title) +
  scale_color_gradientn(
    colors=rev(c('#2b8cbe', 'grey', '#e38071', '#e34a33', '#e31e00')),
    breaks=c(0.05, 0.01, 0.001, 0.0001),
    limits=c(0.00001, 1), trans='log10', oob = scales::squish
  ) +
  scale_y_discrete(name="", limits=rev(topsets), labels=rev(topsets)) +
  theme(
    axis.text.x=element_text(angle=30, hjust=1),
    axis.text.y=element_text(size=8)
  )
}

multi_hitlist_enrichdot <- function(enrichment, plot_title=opt$plot_title){

  #' Create a dotplot showing top enriched gene sets (pathways) for multiple hitlists.
  #' 
  #' @description A ClusterProfiler enrichment result for multiple hitlists is visualized
  #' side-by-side, as dotplot. 
  #' 
  #' @param enrichment. Result of clusterProfiler::enricher.
  #' @param plot_title string. Title of the figure.
  #' @usage single_hitlist_enrichdot(enrichment, plot_title="Top gene sets")
  #' @return ggplot
  #' @details Color corresponds to significance, while size shows gene count in hit list.
  
  #' @examples
  #' single_hitlist_enrichdot(enrichment)
  #' single_hitlist_enrichdot(enrichment, plot_title="Top gene sets")

resulTab <- enrichment@compareClusterResult
topsets <- rownames(resulTab[!duplicated(resulTab$ID), ])[1:30]
hitnames <- levels(resulTab$group)
clusterProfiler::dotplot(enrichment, showCategory=30) +
  labs(title=plot_title) +
  scale_color_gradientn(
    colors=rev(c('#2b8cbe', 'grey', '#e38071', '#e34a33', '#e31e00')),
    breaks=c(0.05, 0.01, 0.001, 0.0001),
    limits=c(0.00001, 1), trans='log10', oob = scales::squish
  ) +
  scale_y_discrete(name="", limits=rev(topsets), labels=rev(topsets)) +
  scale_x_discrete(name="", labels=hitnames) +
  theme(
    axis.text.x=element_text(angle=30, hjust=1),
    axis.text.y=element_text(size=8)
  )
}

single_hitlist_genedot <- function(enrichment){
  
  #' Create a dotplot showing gene set membership of top (best known) genes.
  #' 
  #' @description The top gene sets are scanned for hit genes, until 25 hit genes are
  #' collected. Membership of these top genes is shown in top gene sets.
  #' 
  #' @param enrichment ClusterProfiler result object. Result of an enrichment analysis.
  #' @usage single_hitlist_genedot(enrichment)
  #' @return ggplot
  #' @details Dot size shows how many genes the gene set consists of. Color shows how 
  #' gene sets linked to a certain gene ranked.
 
  #' @examples
  #' single_hitlist_genedot(enrichment)

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
    theme(
      axis.text.x=element_text(size=7, angle=30, hjust=1),
      axis.text.y=element_text(size=8),
      legend.position="bottom",
      legend.justification = c(0,1),
      legend.margin = margin(l=-90, r=90, unit="pt")
    )
}

multi_hitlist_genedot <- function(enrichment, cohort_order=NULL){
  
  #' Create a dotplot showing gene set membership of top (best known) genes.
  #' 
  #' @description The top gene sets of multi-comparison are scanned for hit genes, until
  #' 25 hit genes are collected. Membership of these top genes is shown in top gene sets.
  #' 
  #' @param enrichment ClusterProfiler result object. Result of an enrichment analysis.
  #' @usage multi_hitlist_genedot(enrichment)
  #' @return ggplot
  #' @details Dot size shows how many genes the gene set consists of. Color shows how 
  #' gene sets linked to a certain gene ranked.
 
  #' @examples
  #' multi_hitlist_genedot(enrichment)

  
  enrichment <- enrichment@compareClusterResult
  if(is.null(cohort_order)){
    cohort_order <- unique(as.character(enrichment$group))
  }

  topenr <- enrichment %>%
    .[.$ID %in% .$ID[!duplicated(.$ID)][1:30], ] %>%
    .[, c("ID", "p.adjust", "BgRatio", "geneID")]

  detailedsets <- c()
  detailedgenes <- c()
  for (geneset in unique(topenr$ID)){
    if(length(detailedgenes) < 25){
      setDats <- topenr[topenr$ID == geneset, "geneID"]
      detailedgenes <- unique(c(detailedgenes, unlist(sapply(setDats, strsplit, "/"))))
      detailedsets <- c(detailedsets, geneset)
    }
  }

  geneFuns <- enrichment %>%
    .[detailedsets, c("ID", "p.adjust", "BgRatio", "geneID", "group")] %>%
    transform(BgRatio = sapply(BgRatio, function(x){unlist(strsplit(x, "/"))[[1]]})) %>%
    transform(GeneSetSize = as.numeric(BgRatio)) %>%
    transform(ScaledGeneSetSize = log10(GeneSetSize)) %>%
    transform(ScaledGeneSetSize = ScaledGeneSetSize/min(ScaledGeneSetSize, na.rm=TRUE)) %>%
    transform(group = as.factor(group), levels=cohort_order) %>%
    transform(geneID = as.character(geneID)) %>%
    transform(geneID = strsplit(geneID, "/")) %>%
    transform(geneSet = as.numeric(as.factor(ID))) %>%
    unnest(geneID) %>%
    transform(geneSet = as.numeric(as.factor(ID))) %>%
    transform(gene = as.numeric(as.factor(geneID))) %>%
    transform(cnt = 1) %>%
    pivot_wider(values_from=cnt, names_from=group, values_fill=0)

  ggplot() +
    geom_scatterpie(aes(x=gene, y=geneSet, r=ScaledGeneSetSize), data=geneFuns, cols=cohort_order) +
    scale_x_continuous(breaks=seq(1, length(unique(geneFuns$geneID))), labels = unique(arrange(geneFuns, by=gene)$geneID), "") +
    scale_y_continuous(breaks=seq(1, length(unique(geneFuns$ID))), labels = unique(arrange(geneFuns, by=geneSet)$ID), "")
    geom_scatterpie_legend(geneFuns$GeneSetSize, x=5, y=5) +
    coord_equal() +
    theme(
      axis.text.x=element_text(size=7, angle=30, hjust=1),
      axis.text.y=element_text(size=8),
      legend.position="bottom",
      legend.justification = c(0,1),
      legend.margin = margin(l=-90, r=90, unit="pt")
    )
}

# Ensuring command line connectivity by sourcing an argument parser
source("commandR.r")