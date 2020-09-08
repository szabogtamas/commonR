#!/usr/bin/env Rscript

# A script that takes a hitlist and shows top gene ontologies

library(tidyr)
library(ggplot2)
library(msigdbr)
library(clusterProfiler)

args = commandArgs(trailingOnly=TRUE)

outFile <- args[1]
hitGenes <- unlist(strsplit(args[2], ",", fixed=TRUE))


geneSet <-msigdbr(species="Mus musculus", category="C5", subcategory="BP")
geneSet$gs_name <- gsub('GO_', '', geneSet$gs_name)
geneSet$gs_name <- gsub('_', ' ', geneSet$gs_name)

pdf(outFile)

enricher(hitGenes, TERM2GENE=geneSet[,c('gs_name', 'gene_symbol')], pAdjustMethod="none",)  %>%
  clusterProfiler::dotplot(showCategory=30) + labs(title='Gene ontologies\nassociated with\nhit genes') +
  scale_color_gradientn(colors=rev(c('#2b8cbe', 'grey', '#e38071', '#e34a33', '#e31e00')), breaks=c(0.05, 0.01, 0.001, 0.0001), limits=c(0.00001, 1), trans='log10', oob = scales::squish) +
  theme(axis.text.x=element_text(angle=30, hjust=1))

dev.off()