# A script that plots Kaplan-Meier curves from a data table containing
# both clinical information (survival) and gene expression.

library(optparse)
library(docstring)
library(survival)
library(survminer)

print(commandArgs())

plotKMpair <- function(
  survivaltab,
  ptitle,
  cpalette=c('#e34a33', '#2b8cbe'),
  shortlabels=c("Low", "High"),
  survtype="Overall survival (months)"
  )
#' Plot two Kaplan-Meier curves on a single plot.
#' 
#' @description This plots a Kaplan-Meier for two (sub)populations
#' on a single figure in order to enable comparison of survival
#' under a given condition and in its absence.
#' 
#' @param survivaltab dataframe. Data to plot as a table with "event" and "time"
#' columns, plus a grouping column called "label".
#' @param ptitle character. Title of the KM plot.
#' @param cpalette vector, optional. Colors to be used. Blue and red by default.
#' @param shortlabels vector, optional. Row names for the summary table.
#' @param survtype character. Title of the x axis, telling what kind of survival is shown.
#' @usage plotKMpair(survivaltable, plot_title)
#' @return The figure with two KM curves and a summary table.
#' @details The survival table has to have an event, a time and a label columns.
#' @examples
#' plotKMpair(survivaltable, "A survival example")
#' plotKMpair(survivaltable, cpalette=c('#e34a33', '#2b8cbe'))
#' plotKMpair(survivaltable, cpalette=c('#e34a33', '#2b8cbe'), shortlabels=c("Low", "High"), survtype="Disease-free interval (months)")
#' @export

{
  svive <- Surv(survivaltab$time, survivaltab$event)
  mfit <- do.call(survfit, list(formula=svive ~ label, data=survivaltab))
  survplot <- ggsurvplot(fit=mfit, xlab=survtype, ylab="Survival (%)", surv.scale='percent',
                         legend.title=ptitle, legend=c(0.7, 0.75), tables.col='strata', risk.table=TRUE, pval=TRUE,
                         tables.height=0.2, fontsize=4, palette=cpalette, legend.labs=levels(survivaltab$label))
  survplot$table <- survplot$table +
    theme(axis.text.y = element_text(color="black"), axis.ticks.y = element_blank(), axis.title.y=element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x=element_blank(),
          axis.line=element_blank(), plot.title=element_text(size=12)) +
    scale_y_discrete(labels=shortlabels) +
    theme(axis.text.y = element_text(hjust=0))
  survplot$plot <- survplot$plot + theme(legend.background=element_blank())
  return(survplot)
}

showKMpairs <- function(clinicals, genes, gexprefix="gex_") 
#' Plot two Kaplan-Meier curves on a single plot for a series of genes.
#' 
#' @description This plots a Kaplan-Meier for both the low and the high
#' expressing population for a given gene and arranges these plots into
#' a single figure containing information on a series of genes.
#' 
#' @param dt dataframe. Data to plot as a table with "event" and "time"
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
  n <- 0
  survivalplots <- list()
  for (gene in genes){
    n <- n + 1
    gn <- paste0(gexprefix, gene)
    survivaltab <- clinicals[, c('DFI', 'DFI.time', gn)]
    colnames(survivaltab) <- c('event', 'time', 'gex')
    survivaltab <- survivaltab[complete.cases(survivaltab),]
    mgex <- median(survivaltab$gex)
    survivaltab$label <- ifelse(survivaltab$gex < mgex, 'Low expression', 'High expression')
    cohort <- "BRCA"
    ptitle <- paste0(gene, " expression in ", cohort)
    genesym <- gene #genes_renamed[[gene]]
    survplot <- plotKMpair(survivaltab, ptitle)
    print(survplot)
    survivalplots[[n]] <- survplot
  }
  surv_grob <- arrange_ggsurvplots(survivalplots, nrow=2, ncol=3, print=FALSE)
  p <- ggarrange(plotlist=surv_grob)
  return(p)
}