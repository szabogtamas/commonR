#!/usr/bin/env Rscript

# A script that takes a hitlist and shows top gene ontologies

library(tidyr)
library(ggplot2)
library(msigdbr)
library(clusterProfiler)

args = commandArgs(trailingOnly=TRUE)

outFile <- args[1]
hitGenes <- unlist(strsplit(args[2], ",", fixed=TRUE))

msig_species <- "Mus musculus"
msig_category <- "C5"
msig_subcategory <- "BP"

geneSet <-msigdbr(species=msig_species, category=msig_category, subcategory=msig_subcategory)
geneSet$gs_name <- gsub('GO_', '', geneSet$gs_name)
geneSet$gs_name <- gsub('_', ' ', geneSet$gs_name)

hitGenes  %>%
  enricher(TERM2GENE=geneSet[,c('gs_name', 'gene_symbol')], pAdjustMethod="none", qvalueCutoff=1) %>%
  clusterProfiler::dotplot(showCategory=30) +
  labs(title='Biological processes associated') +
  scale_color_gradientn(
    colors=rev(c('#2b8cbe', 'grey', '#e38071', '#e34a33', '#e31e00')),
    breaks=c(0.05, 0.01, 0.001, 0.0001),
    limits=c(0.00001, 1), trans='log10', oob = scales::squish) +
  theme(axis.text.x=element_text(angle=30, hjust=1))

pdf(outFile)

enricher(hitGenes, TERM2GENE=geneSet[,c('gs_name', 'gene_symbol')], pAdjustMethod="none",)  %>%
  clusterProfiler::dotplot(showCategory=30) + labs(title='Gene ontologies\nassociated with\nhit genes') +
  scale_color_gradientn(colors=rev(c('#2b8cbe', 'grey', '#e38071', '#e34a33', '#e31e00')), breaks=c(0.05, 0.01, 0.001, 0.0001), limits=c(0.00001, 1), trans='log10', oob = scales::squish) +
  theme(axis.text.x=element_text(angle=30, hjust=1))

dev.off()