#!/usr/bin/env Rscript

scriptDescription <- "A script that plots Kaplan-Meier curves from a table containing both clinical information (survival) and gene expression."

scriptMandatoryArgs <- list(
  survivalTab = list(
    abbr="-i",
    type="table",
    readoptions=list(sep="\t", stringsAsFactors=FALSE),
    help="Table of survival times, event and gene expression."
  )
)

scriptOptionalArgs <- list(
  eventLabel = list(
    default="OS",
    abbr="-e",
    help="Column name for event information."
  ),
  timeLabel = list(
    default="OS.time",
    abbr="-t",
    help="Column name for survival time information."
  ),
  gexLabel = list(
    default="gex",
    abbr="-x",
    help="Column name for gene expression."
  ),
  conditionColors = list(
    default=NULL,
    type="vector",
    help="Colors for points belonging to a given condition on PCA."
  ),
  commandRpath = list(
    default="commandR.r",
    help="Path to command line connectivity script (if not in cwd)."
  )
)

opt <- list()
for (rn in names(scriptOptionalArgs)){
  opt[[rn]] <- scriptOptionalArgs[[rn]][["default"]]
}

for (pk in c("tidyr", "dplyr", "purrr", "tibble", "survival", "survminer")){
  if(!(pk %in% (.packages()))){
    library(pk, character.only=TRUE)
  }
}

#' The main function of the script, executed only if called from command line.
#' Calls subfunctions according to supplied command line arguments.
#' 
#' @param opt list. a named list of all command line options; will be passed on 
#' 
#' @return Not intended to return anything, but rather save outputs to files.
main <- function(opt){
  
  opt$help <- NULL
  opt$verbose <- NULL

  cat("Plotting Kaplam-Meier curves\n")
  do.call(plotKMpair, opt)

  invisible(NULL)
}


plotKMpair <- function(
  survivaltab,
  ptitle,
  cpalette=c('#e34a33', '#2b8cbe'),
  shortlabels=c("Low", "High"),
  survtype="Overall survival (months)"
)
{
  #' Plot two Kaplan-Meier curves on a single plot.
  #' 
  #' @description This plots Kaplan-Meier curves for (sub)populations
  #' on a single figure in order to enable comparison of survival
  #' under a given condition and in its absence.
  #' 
  #' @param survivaltab dataframe. Data to plot as a table with "event" and "time"
  #' columns, plus a grouping column called "label".
  #' @param ptitle character. Title of the KM plot.
  #' @param cpalette vector, optional. Colors to be used. Blue and red by default.
  #' @param shortlabels vector, optional. Row names for the summary table.
  #' @param survtype character, optional. Title of the x axis, telling what kind
  #' of survival is shown.
  #' @usage plotKMpair(survivaltable, plot_title)
  #' @return The figure with two KM curves and a summary table.
  #' @details The survival table has to have event, time and label columns.
  #' By default, it calculates overall survial.
  #' @examples
  #' survivaltable <- data.frame(event=c(NA, NA, 0, 1, 1, 1), time=c(NA, NA, 1556, 490, 3014, 2165), gex=c(17.63, 16.17, 17.79, 18.26, 17.04, 17.69), label=c('Low expression', 'Low expression', 'High expression', 'High expression', 'High expression', 'Low expression'))
  #' plotKMpair(survivaltable, "A survival example")
  #' plotKMpair(survivaltable, cpalette=c('#e34a33', '#2b8cbe'))
  #' plotKMpair(survivaltable, cpalette=c('#e34a33', '#2b8cbe'), shortlabels=c("Low", "High"), survtype="Disease-free interval (months)")
  
  svive <- Surv(survivaltab$time, survivaltab$event)
  mfit <- do.call(survfit, list(formula=svive ~ label, data=survivaltab))
  survplot <- ggsurvplot(fit=mfit, xlab=survtype, ylab="Survival (%)", surv.scale='percent',
                         legend.title=ptitle, legend=c(0.6, 0.4), tables.col='strata', risk.table=TRUE, pval=TRUE,
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

showKMpairs <- function(
  clinicals,
  genes,
  gexprefix='gex_',
  eventcol='DFI',
  timecol='DFI.time',
  timefactor=30,
  cohort=NULL,
  title_text=NULL,
  genedict=NULL,
  numrows=2,
  numcols=3,
  cpalette=c('#e34a33', '#2b8cbe'),
  shortlabels=c('Low', 'High'),
  survtype='Overall survival (months)'
)
{
  #' Plot two Kaplan-Meier curves on a single plot for a series of genes.
  #' 
  #' @description This plots a Kaplan-Meier for both the low and the high
  #' expressing population for a given gene and arranges these plots into
  #' a single figure containing information on a series of genes.
  #' 
  #' @param cliniclas dataframe. Clinicla data contaning survival information,
  #' as well as gene expression values.
  #' @param genes vector. Genes to be shown on the figure.
  #' @param gexprefix character, optional. Text to be appended to the gene name
  #' and used as title.
  #' @param eventcol character, optional. Column name for the event.
  #' @param timecol character, optional. Column name for time untill the event.
  #' @param timefactor numeric, optional. Factor to scale time data.
  #' Default is 30, converting days to months.
  #' @param cohort character, optional. Cohort or diagnosis to be added to the title.
  #' @param title_text character, optional. Text to be appended to the gene name
  #' and used as title.
  #' @param genedict list, optional. Gene names be used intead of symbols.
  #' @param numrows numeric, optional. Number of subfigures in a row.
  #' @param numcols numeric, optional. Number of subfigures in a column.
  #' @param cpalette vector, optional. Colors to be used. Blue and red by default.
  #' @param shortlabels vector, optional. Row names for the summary table.
  #' @param survtype character, optional. Title of the x axis, telling what kind
  #' of survival is shown.
  #' @usage showKMpairs(clinicaltable, genes)
  #' @return The figure with two KM curves and a summary table.
  #' @details The survival table has to have an event, a time and a label columns.
  #' By default, it calculates overall survial.
  #' @examples
  #' clinicaltable <- data.frame(sample=c("TCGA-BH-A1ES-06A", "TCGA-BH-A1FE-06A", "TCGA-BH-A18V-06A", "TCGA-B6-A0X1-01A", "TCGA-GM-A2DA-01A", "TCGA-B6-A0RH-01A"), type=c("BRCA", "BRCA", "BRCA", "BRCA", "BRCA", "BRCA"), DFI.time=c(NA, NA, 1556, 490, 3014, 2165), DFI=c(NA, NA, 0, 1, 1, 1), gex_TP53=c(17.63, 16.17, 17.79, 18.26, 17.04, 17.69), gex_PTEN=c(9.903, 9.867, 17.110, 16.080, 14.380, 16.630))
  #' showKMpairs(clinicaltable, c('TP52', 'PTEN'))
  #' showKMpairs(clinicaltable, c('TP52', 'PTEN'), gexprefix='gex_', eventcol='DFI', timecol='DFI.time', timefactor=30, cohort='BRCA', title_text=" expressed", genedict=list(TP53='p53'), numrows=2, numcols=3, cpalette=c('#e34a33', '#2b8cbe'), shortlabels=c("Low", "High"), survtype="Disease-free interval (months)")
  
  n <- 0
  survivalplots <- list()
  for (gene in genes){
    n <- n + 1
    gn <- paste0(gexprefix, gene)
    survivaltab <- clinicals[, c(eventcol, timecol, gn)]
    colnames(survivaltab) <- c('event', 'time', 'gex')
    survivaltab$time <- survivaltab$time/timefactor
    survivaltab <- survivaltab[complete.cases(survivaltab),]
    mgex <- median(survivaltab$gex)
    survivaltab$label <- ifelse(survivaltab$gex < mgex, 'Low expression', 'High expression')
    
    genesym <- gene
    if(is.null(genesym)){
      if(gene %in% names(genedict)){
        genesym <- genedict[[gene]]
      }
    }
    if(is.null(title_text)){
      title_text <- " expression"
    }
    if(is.null(cohort)){
      ptitle <- paste0(genesym, title_text)
    } else {
      ptitle <- paste0(genesym, title_text, " in ", cohort)
    }
    
    survplot <- plotKMpair(survivaltab, ptitle, cpalette=cpalette, shortlabels=shortlabels, survtype=survtype)
    survivalplots[[n]] <- survplot
  }
  surv_grob <- arrange_ggsurvplots(survivalplots, nrow=numrows, ncol=numcols, print=FALSE)
  return(ggarrange(plotlist=surv_grob))
}

# Ensuring command line connectivity by sourcing an argument parser
source(opt$commandRpath, local=TRUE)