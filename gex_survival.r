# A script that plots Kaplan-Meier curves from a data table containing
# both clinical information (survival) and gene expression

library(optigrab)
library(survminer)

opts <- c( "--fn", "TCGA-BRCA_data.tsv")
opt_get( "fn", opts=opts)
opt_help()

clinicals <- read.table(fn, header = FALSE, sep = "\t")
print(head(clinicals))
