# A script that plots Kaplan-Meier curves from a data table containing
# both clinical information (survival) and gene expression.

library(optparse)
library(docstring)
library(survival)
library(survminer)

# Parse command line options
parser <- OptionParser()
parser <- add_option(parser, c("-v", "--verbose"), action="store_true", default=FALSE,
                     help="Print some progress messages to stdout.")
parser <- add_option(parser, c("-q", "--quietly"), action="store_false", 
                     dest="verbose", help="Create figures quietly, without printing to stdout.")
parser <- add_option(parser, c("-f", "--inputfile"), default=NULL, 
                     help="Data file containing survival information and gene expression.",
                     metavar="path_to_clinicals")
parser <- add_option(parser, c("-g", "--genes"), default=NULL, 
                     help="Data file containing survival information and gene expression.",
                     metavar="genelist")
parser <- add_option(parser, c("--km_out"), default="survival.pdf", 
                     help="Path to save the Kaplan-Meier plots default: [%default]",
                     metavar="kaplanmeier_file")
parser <- add_option(parser, c("--event"), default="OS", 
                     help="Column name for the event: [%default]",
                     metavar="survival event column")
parser <- add_option(parser, c("--time"), default="OS.time", 
                     help="Column for survival times: [%default]",
                     metavar="survival time column")
parser <- add_option(parser, c("--timefactor"), default=30, type="integer", 
                     help="Scaling factor for time data (e.g. converting to months): [%default]",
                     metavar="time factor")
parser <- add_option(parser, c("--cohort"), default=NULL, 
                     help="Cohort name: [%default]")
parser <- add_option(parser, c("--title_text"), default=NULL, 
                     help="Text added to subplot titles: [%default]",
                     metavar="title text")
parser <- add_option(parser, c("--genedict"), default=NULL, 
                     help="A gene name mapping: [%default]",
                     metavar="gene_names")
parser <- add_option(parser, c("--numrows"), default=2, type="integer", 
                     help="Number of subplots in a row: [%default]")
parser <- add_option(parser, c("--numcols"), default=3, type="integer", 
                     help="Number of subplots in a column: [%default]")
parser <- add_option(parser, c("--cpalette"), default="#e34a33,#2b8cbe", 
                     help="Color palette: [%default]",
                     metavar="colors")
parser <- add_option(parser, c("--shortlabels"), default="Low,High", 
                     help="Short labels that with appear in the survival table: [%default]",
                     metavar="table_rows")
parser <- add_option(parser, c("--survtype"), default="Disease-free interval (months)", 
                     help="Path to save the Kaplan-Meier plots default: [%default]",
                     metavar="x_label")
parser <- add_option(parser, "--gexprefix", default="gex_", metavar="column name prefix",
                     help="Prefix for the gene expression column: [%default]")
opt <- parse_args(parser)


# Check if mandatory arguments are present
if ( is.null(opt$inputfile) ) { 
  if ( opt$verbose ) { 
    write("Sorry, cannot proceed without a data table. Please provide a path to survival and gene expression data...\n", stderr())
  }
  checkpass <- FALSE
} else {
  checkpass <- TRUE
}
if ( is.null(opt$genes) ) {  
  if ( opt$verbose ) { 
    write("No genes given. Cannot guess what genes are to be analysed.\n", stderr())
  }
  checkpass <- FALSE
} else {
  opt$genes <- unlist(strsplit(opt$genes, ","))
  checkpass <- TRUE
}
if ( is.null(opt$genedict) ) {  
  opt$genedict <- list()
} else {
  ntab <- read.table(opt$genedict, header=FALSE, sep="\t", row.names=2)
  opt$genedict <- ntab[,1]
  names(opt$genedict) <- rownames(ntab)
}

opt$cpalette <- unlist(strsplit(opt$cpalette, ","))
opt$shortlabels <- unlist(strsplit(opt$shortlabels, ","))

# Define helper functions

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

# Execute main function if mandatory arguments are supplied via the command line (otherwise print help message)
if ( checkpass ) { 
  clinicals <- read.table(opt$inputfile, header=TRUE, sep="\t")
  pdf(file=opt$km_out,7.2, 5.4, onefile=FALSE)
  print(showKMpairs(clinicals, opt$genes, gexprefix=opt$gexprefix, eventcol=opt$eventcol, timecol=opt$timecol, timefactor=opt$timefactor, cohort=opt$cohort, title_text=opt$title_text, genedict=opt$genedict, numrows=opt$numrows, numcols=opt$numcols, cpalette=opt$cpalette, shortlabels=opt$shortlabels, survtype=opt$survtype))
  if ( opt$verbose ) { 
    cat(paste0("Kaplan-Meier plots saved as ", opt$km_out, "\n")) 
  }
  dev.off()
} else {
  print_help(parser)
}