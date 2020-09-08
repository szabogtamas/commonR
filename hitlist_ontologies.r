#!/usr/bin/env Rscript

# A script that takes a hitlist and shows top gene ontologies

library(tidyr)
library(ggplot2)
library(msigdbr)
library(clusterProfiler)

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
  p2 <- cnetplot(enrichment)

  grid.arrange(p1, p2, ncol=1)

  pdf(outFile, height=9.6, width=7.2)
  print(p1)
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
    limits=c(0.00001, 1), trans='log10', oob = scales::squish) +
  theme(axis.text.x=element_text(angle=30, hjust=1))
}

if (!interactive()) {
  main()
}