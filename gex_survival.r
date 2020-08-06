# A script that plots Kaplan-Meier curves from a data table containing
# both clinical information (survival) and gene expression.

library(docstring)
library(optigrab)
library(survival)
library(survminer)

fn <- opt_get( "fn", default="TCGA-BRCA_data.tsv")
cohort <- opt_get( "cohort", default="")
genes <- opt_get( "genes", default="TP53")
gendict <- opt_get( "gendict", default=NULL)
timecol <- opt_get( "timecol", default="OS.time")
eventcol <- opt_get( "eventcol", default="OS")
gexprefix <- opt_get( "gexprefix", default="gex_")
col2plot <- opt_get( "col2plot", default="3")
row2plot <- opt_get( "row2plot", default="2")

plotKMpair <- function(dt) 
#' Plot two Kaplan-Meier curves on a single plot.
#' 
#' @description This plots a Kaplan-Meier for two (sub)populations
#' on a single figure in order to enable comparison of survival
#' under a given condition and in its absence.
#' 
#' @param dt dataframe. Data to plot as a tble with "event" and "time"
#' columns, plus a grouping column
#' @param y character. The second item to paste Defaults to "!" but
#' "?" would be pretty great too
#' @usage mypaste(x, y)
#' @return The figure with two KM curves and a summary table.
#' @details The inputs can be anything that can be input into
#' the paste function.
#' @examples
#' mypaste(1, 3)
#' mypaste("hey", "you")
#' mypaste("single param")
#' @export
#' @importFrom base survminer
{
  dualcolors <- c('#e34a33', '#2b8cbe')
  svive <- Surv(survivaltab$time, survivaltab$event)
  mfit <- do.call(survfit, list(formula=svive ~ label, data=survivaltab))
  #ptitle <- paste(c(genes_renamed[[gene]], " expression"), collapse = '')
  ptitle <- 'Halo'
  msurvplot <- ggsurvplot(fit=mfit, xlab="Disease-free interval (months)", ylab="Survival (%)", surv.scale='percent', legend.title = ptitle, legend=c(0.7, 0.75), risk.table=TRUE, tables.height=0.2, fontsize=4, palette=dualcolors, pval=TRUE, legend.labs=levels(survivaltab$label), tables.col='strata')
  msurvplot$table <- msurvplot$table + theme(axis.text.y = element_text(color="black"), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.line=element_blank(), plot.title=element_text(size=12))
  #msurvplot$table <- msurvplot$table + scale_y_discrete(labels=c("Low", "High")) + theme(axis.text.y = element_text(hjust=0))
  msurvplot$table <- msurvplot$table + scale_y_discrete(labels=c('', '')) + theme(legend.background=element_blank())
}

showKMpairs <- function(fn, genes, gexprefix="gex_") 
#' Plot two Kaplan-Meier curves on a single plot for a series of genes.
#' 
#' @description This plots a Kaplan-Meier for both the low and the high
#' expressing population for a given gene and arranges these plots into
#' a single figure containing information on a series of genes.
#' 
#' @param dt dataframe. Data to plot as a tble with "event" and "time"
#' columns, plus a grouping column
#' @param y character. The second item to paste Defaults to "!" but
#' "?" would be pretty great too
#' @usage mypaste(x, y)
#' @return The figure with two KM curves and a summary table.
#' @details The inputs can be anything that can be input into
#' the paste function.
#' @examples
#' mypaste(1, 3)
#' mypaste("hey", "you")
#' mypaste("single param")
#' @export
#' @importFrom base survminer
{
  clinicals <- read.table(fn, header=TRUE, sep="\t")
  n <- 0
  survivalplots <- list()
  for (gene in genes){
    n <- n + 1
    gn <- paste0(gexprefix, gene)
    survivaltab <- clinicals[, c('DFI', 'DFI.time', gn)]
    colnames(survivaltab) <- c('event', 'time', 'gex')
    mgex <- median(survivaltab$gex)
    survivaltab$label <- ifelse(survivaltab$gex < mgex, 'Low expression', 'High expression')
    plotKMpair(survivaltab)
    survivalplots[[n]] <- msurvplot
    print(msurvplot)
  }
  surv_grob <- arrange_ggsurvplots(survivalplots, nrow=2, ncol=3, print=FALSE)
  p <- ggarrange(plotlist=surv_grob)
  return(p)
}

pdf(file=outname)
showKMpairs(fn, genes=c("TP53", "PIDD1"))
dev.off()
print("hey")
