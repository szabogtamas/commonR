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
  outFile = list(
    abbr="-o",
    default=NULL,
    help="Prefix for output files."
  ),
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
  legendLabel = list(
    default="Gene expression",
    help="Legend label for groups."
  ),
  survivalLabel = list(
    default="Disease-free interval (months)",
    help="Label for x axis (survival type and time)."
  ),
  conditionLabels = list(
    default=NULL,
    type="vector",
    help="Labels for the two conditions (high and low expression)."
  ),
  shortLabels = list(
    default=NULL,
    type="vector",
    help="Short version of condition labels for risk table."
  ),
  conditionColors = list(
    default=NULL,
    type="vector",
    help="Colors for conditions."
  ),
  pctHigh = list(
    default=0.75,
    type="vector",
    help="Percentile limit for high expression."
  ),
  pctLow = list(
    default=0.25,
    type="vector",
    help="Percentile limit for low expression."
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

for (pk in c("tidyr", "dplyr", "survival", "survminer")){
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

  outFile <- opt$outFile
  if(is.null(outFile)){
    outFile <- paste(opt$gexLabel, "_survival")
  }
  opt$outFile <- NULL
  opt$help <- NULL
  opt$verbose <- NULL

  cat("Plotting Kaplan-Meier curves\n")
  p <- do.call(plotKMpair, opt)

  cat("Saving figure\n")
  fig2pdf(p, paste0(outFile, ".pdf"), height=3.6, width=7.2)

  invisible(NULL)
}

#' Plot two Kaplan-Meier curves on a single plot.
#' 
#' @description This plots two Kaplan-Meier curves for (sub)populations
#' on a single figure in order to enable comparison of survival
#' under two conditions (typically high and low expression of a gene).
#' 
#' @param survivalTab dataframe. Data to plot as a table with "event" and "time"
#' columns, plus gene expression information.
#' @param eventLabel character. Column name for event information.
#' @param timeLabel character. Column name for time information.
#' @param gexLabel character, optional. Row names for the summary table.
#' @param legendLabel character, optional. Column name for gene expression.
#' @param survivalLabel character, optional. Label for x axis (survival type and time).
#' @param conditionLabels vector, optional. Labels for the two conditions (high and low expression).
#' @param shortLabels vector, optional. Short version of condition labels for risk table.
#' @param conditionColors vector, optional. Colors for conditions.
#' @param pctHigh float, optional. Percentile limit for high expression.
#' @param pctLow float, optional. Percentile limit for low expression.
#' @param ... ellipse, optional. Arguments passed additionally to ggsurvplot.
#' @usage plotKMpair(survivalTab, eventLabel, timeLabel)
#' @return The figure with two KM curves and a summary table.
#' 
#' @examples
#' plotKMpair(survivalTab, eventLabel, timeLabel, conditionColors=c('#e34a33', '#2b8cbe'), shortlabels=c("Low", "High"), survivalLabel="Disease-free interval (months)")
plotKMpair <- function(
  survivalTab,
  eventLabel,
  timeLabel,
  gexLabel = "gex",
  legendLabel = "Gene expression",
  survivalLabel = "Survival in months",
  conditionLabels = c("High", "Low"),
  shortLabels = NULL,
  conditionColors = NULL,
  pctHigh= 0.5,
  pctLow = 0.5,
  ...
){
  
  if (is.null(conditionColors)){
      conditionColors <- default_colors
    }
  
  if(is.null(shortLabels)){
    shortLabels <- conditionLabels
  }

  ### Subset input table for columns to be used and standardize column names
  stdSurvival <- survivalTab %>%
    select(one_of(eventLabel, timeLabel, gexLabel)) %>%
    rename(setNames(c(eventLabel, timeLabel, gexLabel), c("event", "time", "gex"))) %>%
    mutate(
      gex_highpct = quantile(gex, pctHigh),
      gex_lowpct = quantile(gex, pctLow),
      label = case_when(
        gex > gex_highpct ~ conditionLabels[1],
        gex <= gex_lowpct ~ conditionLabels[2],
        TRUE ~ NA_character_
      ),
      time = as.numeric(time),
      survival = Surv(time, event)
    ) %>%
    filter(!is.na(label))

  ### Plot survival
  survplot <- survfit(formula=survival ~ label, data=stdSurvival) %>%
    ggsurvplot(
      xlab=survivalLabel,
      ylab="Survival (%)",
      surv.scale='percent',
      legend.title = legendLabel,
      legend=c(0.7, 0.75),
      risk.table=TRUE,
      tables.height=0.2,
      fontsize=4,
      palette=conditionColors,
      legend.labs=conditionLabels,
      pval=TRUE,
      tables.col='strata',
      ...
    )
  
  ### Optimize appearance of plot
  survplot$table <- survplot$table +
    scale_y_discrete(labels=shortLabels) +
    theme(
      axis.text.y = element_text(color="black", hjust=0),
      axis.ticks.y = element_blank(),
      axis.title.y=element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x=element_blank(),
      axis.line=element_blank(),
      plot.title=element_text(size=12)
    )


  survplot$plot <- survplot$plot +
    theme(legend.background=element_blank())
  
  return(survplot)
}

# Ensuring command line connectivity by sourcing an argument parser
source(opt$commandRpath, local=TRUE)