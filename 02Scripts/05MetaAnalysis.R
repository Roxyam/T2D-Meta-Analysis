#!/usr/bin/env Rscript

# ~~~~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~ ** ~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~~~~ #
#
# 05MetaAnalysis.R
#
#
# Computes the meta-analysis for a given set of studies:
#
#     This script computes the meta-analysis of a study for a given 
#     contrast using a random effects model and the DL method.
#     It also generates a report with the results.
#
#
#     ****************************************************************
#     *    Author:  Roxana Andreea Moldovan                          *
#     *    Contact: roxana.andreea.moldovan@gmail.com                *
#     *    Version: 1.0                                              *
#     *    Creation date: 02/04/2023                                 *
#     ****************************************************************
#
#
# ~~~~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~ ** ~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~~~~ #


# ~~~~~~~~~~~~ Loading packages ~~~~~~~~~~~~ #

library(pacman)
pacman::p_load(argparse)
pacman::p_load(glue)
pacman::p_load(limma)
pacman::p_load(metafor)
pacman::p_load(parallel)
pacman::p_load(stringr)


# ~~~~~~~~~~~~ Functions ~~~~~~~~~~~~ #

source("Functions/Functions.R")


# ~~~~~~~~~~~~ Parameters ~~~~~~~~~~~~ #

parser <- ArgumentParser(description="This script compute meta analysis
                                      for a set of studies")
# Studies
parser$add_argument("-s", "--studies",
                    action="store",
                    type="character",
                    default="ALL",
                    help="Studies to include
                          separated by comma,
                          by default all will
                          be taken.")

# Differential expression file
parser$add_argument("-f", "--filein",
                    action="store",
                    type="character",
                    required=TRUE,
                    help="The path and filename for the file
                          that includes all the information
                          of the differential expression analysis
                          to be udes.")

# Tissue
parser$add_argument("-t", "--tissue",
                    action="store",
                    type="character",
                    required=TRUE,
                    help="Tissue or tissues to include.")

# Contrasts
parser$add_argument("-c", "--contrast",
                    action="store",
                    type="character",
                    required=TRUE,
                    help="Contrast to include.")

# Out dir
parser$add_argument("-o", "--outdir",
                    action="store",
                    type="character",
                    default=".",
                    help="Where you would like the output files to be placed,
                          by default the current directory will be taken.")
# Id
parser$add_argument("-id", "--id",
                    action="store",
                    type="character",
                    default="Meta-analysis",
                    help="Job id. By default is 'Meta-analysis'.")
# Output files prefix
parser$add_argument("-p", "--prefix",
                    action="store",
                    type="character",
                    default="Meta-analysis",
                    help="Prefix to use in output files, by default is
                          'Meta-analysis'.")
# Method
parser$add_argument("-m", "--method",
                    action="store",
                    type="character",
                    default="coefSE",
                    choices = c("SE", "coefSE"), 
                    help="Method for SS calculating.")
# Report
parser$add_argument("-r", "--report",
                    action="store",
                    type="character",
                    default=TRUE,
                    help="Create an R Markdown and HTML with the results.")
# Redo
parser$add_argument("-re", "--redo",
                    action="store_true",
                    default = FALSE,
                    help="Repeat the analysis, if it is not specified,
                          the files of previous executions are taken")
# Pvalue
parser$add_argument("-pv", "--plim",
                    action="store",
                    type="integer",
                    default = 0.05,
                    help="Adjusted p-value limit. By default is 0.05.")
# Max
parser$add_argument("-mx", "--pmax",
                    action="store",
                    type="integer",
                    default = 50,
                    help="Maximum number of generated plots.
                    By default is 50")

# ~~~~~~~~~~~~ Main ~~~~~~~~~~~~ #

#------------- Execution 
#args <- parser$parse_args(args = c('-t="SAT"',
#                                   '-f=/home/rmoldovan/T2D-Meta-Analysis/Data/04DE/DifferentialExpressionGroup.RData',
#                                   '-o=/home/rmoldovan/T2D-Meta-Analysis/Data/05MA/Obesity/SDIO',
#                                   '-c=(Ob_IS.M - Np_IS.M) - (Ob_IS.F - Np_IS.F)',
#                                   '-p=Meta-analysis_4',
#                                   '-id=Meta-analysis_4'))
#
#------------- Checking arguments 
# Load the data
Datas = get(load(args$filein))

# Get the studies
if (args$studies != "ALL"){
  # Take studies to include
  Studies = eval(parse(text = args$studies))
  Studies = trimws(Studies, whitespace = "[ \t\r\n\\.]")
} else{
  Studies = names(Datas)
}

# Check
if (! all(Studies %in% names(Datas))){
  cat("The indicated file does not contain information on",
      "the following studies. \n")
  cat("\t", Studies[! Studies %in% names(Datas)])
  cat("\n")
  stop("Please review input.",call. = FALSE)
}

# Get the tissues
Tissue = eval(parse(text = args$tissue))
Tissue = trimws(Tissue, whitespace = "[ \t\r\n\\.]")

# Get the contrasts
Contrasts = stringr::str_split_1(args$contrast, pattern = ",")
Contrasts = trimws(Contrasts, whitespace = "[ \t\r\n\\.]")

# Parameters
method = "DL"
PlotsDir = glue("{args$outdir}/Plots")

# Output directory
if(! file.exists(args$outdir)){
  system(glue("mkdir {args$outdir}"))
}
# Plots directory
system(glue("mkdir {PlotsDir}"))


#------------- 1. Meta-analysis for genes
cat("\n---------------------------------------------\n")
cat("Analyzing...\n") 
contrast = Contrasts
Mats_out = glue("{args$outdir}/{args$prefix}_mats.RData")
Mats_in = glue("{args$outdir}/{args$id}_mats.RData")
MA_out = glue("{args$outdir}/{args$prefix}_MA.RData")
MA_in = glue("{args$outdir}/{args$id}_MA.RData")
DF_out = glue("{args$outdir}/{args$prefix}_DF.RData")
DF_in = glue("{args$outdir}/{args$prefix}_DF.RData")

if(file.exists(Mats_in) & isFALSE(args$redo)){
  # Load previous data
  load(file=Mats_in)
}else{
  Mats = RX_MetaAnalysis_prep(Studies = Studies,
                              Tissue = Tissue,
                              Datas = Datas,
                              SE_coef = "SE",
                              Contrast = contrast)
  save(Mats, file = Mats_out)
}
# Take only present studies
Studies = colnames(Mats$MatLogFC)

if(file.exists(MA_in) & isFALSE(args$redo)){
  # Load previous data
  load(file=MA_in)
}else{
  #------------- 1. Meta-analysis for genes  
  MA <- sapply(rownames(Mats$MatLogFC),
               function(x) metafor::rma(yi = Mats$MatLogFC[x, ],
                                        sei = Mats$MatSE[x, ],
                                        method = method,
                                        slab = Studies,
                                        verbose = FALSE),
               simplify = FALSE,USE.NAMES = TRUE ) 
  save(MA, file = MA_out)
}

if(file.exists(DF_in) & isFALSE(args$redo)){
  # Load previous data
  load(file=DF_in)
}else{
  # Data frame including the detailed results.
  resultMA <- as.data.frame(do.call("rbind",
                                    sapply(MA,
                                           function(x){c("lb" = x$ci.lb, 
                                                         "logFC" = x$b, 
                                                         "ub" = x$ci.ub,
                                                         "pvalue" = x$pval,
                                                         "QE" = x$QE,
                                                         "QEp" = x$QEp, 
                                                         "SE" = x$se,
                                                         "tau2" = x$tau2, 
                                                         "I2" = x$I2, 
                                                         "H2" = x$H2)},
                                           simplify = FALSE)))
  # Adjust p.values
  resultMA$p.adjust.BH <- stats::p.adjust(resultMA$pvalue, method = "fdr")
  resultMA$p.adjust.BY  <- stats::p.adjust(resultMA$pvalue, method = "BY")
  resultMA$logFDR <- -log10(resultMA$p.adjust.BH)
  resultMA$N <- rowSums(!is.na(Mats$MatSE))
  # Add info
  #Significative
  resultMA$Sig = FALSE
  resultMA[which(resultMA$p.adjust.BH < args$plim), "Sig"] = TRUE
  # Influence
  resultMA$N.concordantFC = "-"
  resultMA$incluence.St = "-"
  resultMA$sensi.global = "-"
  resultMA$sensi.specific = "-"
  # Add info
  for (i in which(resultMA$Sig == TRUE)){
    # Studies that include information on this gene
    St <- stringr::str_split_i(colnames(Mats$MatSE)[!is.na(Mats$MatSE)[i,]], pattern = "\\.", i=1) 
    
    
    # Influence info:
    ## Number of studies where the sign of the logFC is the same as the global logFC
    resultMA[i, "N.concordantFC"] <- sum(sign(MA[[i]]$yi) == rep(sign(MA[[i]]$b),length(St)))
    # Influence studies
    inf <- influence(MA[[i]])
    res <- paste(St[inf$is.infl], collapse = ",")
    resultMA[i, "incluence.St"] <- ifelse(res =="", "Non", res)
    
    # Sensitivity info:
    ## Leave one out
    l1 <- as.data.frame(leave1out(MA[[i]]))
    ## Are there differences between global analysis and leaving one study out?
    resultMA[i, "sensi.global"] <- t.test(x= l1$estimate,
                                          mu=as.numeric(MA[[i]]$b))$p.value
    
    ## Changes after leaving one study out
    res2 <- paste(St[l1$pval > 0.05], collapse = ",")
    resultMA[i, "sensi.specific"] <- ifelse(res2 =="", "all.p.values < 0.05", res2)
  }
  # Order
  ord = c("lb", "logFC", "ub", "pvalue", "p.adjust.BH", "p.adjust.BY", "Sig", "N")
  resultMA = data.frame(Symbol = Mats$Anot[rownames(resultMA), "SYMBOL"],
                           resultMA[,ord],
                           resultMA[, which(! colnames(resultMA) %in% ord)]
                           )
  # Save
  save(resultMA, file = DF_out)
}


sigData0 = resultMA[resultMA[, "p.adjust.BH"] < args$plim,]
sigData = sigData0[order(abs(sigData0[, "p.adjust.BH"])),]
if (nrow(sigData) > args$pmax ){
  sigData = sigData[c(1:args$pmax),]
}

color = "#2494b5"
  for (gene in rownames(sigData)){
    res <- MA[[gene]]
    #res <- rma(yi= Mats$MatLogFC[gene,], sei =Mats$MatSE[gene,], method = "DL", slab = Studies)
    # Forest plot
    svg(glue("{PlotsDir}/{args$prefix}_{gene}forest.svg"))
    forest(res, 
           slab = toupper(colnames(Mats$MatLogFC)),
           xlab="Fold Change", cex=0.7,
           mlab="Modelo 'DL' para todos los estudios.", col = color, 
           main = paste("\n", Mats$Anot[gene, "SYMBOL"] , " (", gene, ")",sep=""))    
    text( 9,-3, "logFC [IC 95%]", pos=2, cex = 0.7)
    dev.off()
    
    #Funnel plot
    svg(glue("{PlotsDir}/{args$prefix}_{gene}funnel.svg"))
    funnel(res, main=paste("\n", Mats$Anot[gene, "SYMBOL"] , " (", gene, ")",sep=""), back ="transparent", ylab = "Error estándar",
           xlab = "logFC", label = 2, shade=c("#F0EFEF", "#90C432", "#2494b5"), refline=0,
           level=c(90, 95, 99), legend=TRUE, atransf=exp)
    dev.off()  
    
    #Influence plot
    svg(glue("{PlotsDir}/{args$prefix}_{gene}influence.svg"))
    inf <- influence(res)
    plot(inf)
    dev.off()  
  }

#------------- Report

report_out = glue("{args$outdir}/{args$prefix}Report.rmd")
if(isTRUE(args$report)){
  # File info
  cat('---
title: "05 Metaanálisis."
author: "Roxana Andreea Moldovan Moldovan"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: false
    toc_float: false
    number_sections: false 
    code_folding: hide
    theme: flatly
    bg: "#FFFFFF"
    fg: "#202020"
    primary: "#90C432"
    secondary: "#2494b5" 
  pDT_document:
    extra_dependencies: "subfig"
linkcolor: blue
urlcolor: blue
citecolor: blue
header-includes:
- \\usepackage{subfig}
- \\usepackage{float}
fig_caption: yes 
link-citations: TRUE 
---
        ',
      sep ="",
      file = report_out,
      append = FALSE
  )
  
  ## Meta analisis
  cat('\n##   {.tabset .tabset-fade -} \n',
      '&nbsp;  \n\n',
      sep ="",
      file = report_out,
      append = TRUE)
  
  # Load packages
  cat(
    # Open chunk
    '\n```{r echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE}\n',
    '# Load packages \n',
    'library(pacman)
pacman::p_load(glue)
pacman::p_load(DT)
pacman::p_load(magick)',
    # Close chunk
    '\n```  \n\n',
    sep ="",
    file = report_out,
    append = TRUE)
  
  cat('\nEn este documento se muestran los resultados del metaanálisis de los estudios: ', 
      paste(Studies, collapse = ", "),
      ' para el contraste ',
      args$contrast, 
      '.  \n&nbsp;  \n\n',
      sep ="",
      file = report_out,
      append = TRUE)
  
  cat('\nTomando un p valor ajustado por BH < ', args$plim, ' encontramos ', nrow(sigData0), ' genes significativos.  \n',
      '&nbsp;  \n\n',
      sep ="",
      file = report_out,
      append = TRUE)
  
  cat('\n### Tabla de resultados {-} \n',
      '&nbsp;  \n\n',
      sep ="",
      file = report_out,
      append = TRUE)
  
  cat(
    # Open chunk
    '\n```{r echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE}\n',
    'DT = get(load("', DF_out, '"))\n',
    'datatable(DT, options = list(caption = "Resultados del metaanálisis.",
                              scrollX = TRUE),
                              filter = "top")',
    # Close chunk
    '\n```  \n\n','\n&nbsp;  \n\n',
    sep ="",
    file = report_out,
    append = TRUE)
  
  cat('\n### Representación gráfica {-}  \n&nbsp; \n',
      '\nAquí se muestran las figuras características del metaanálisis para los
genes significativos tomando un p valor ajustado por BH < ', args$plim, '.  \n',
      '&nbsp;  \n\n',
      sep ="",
      file = report_out,
      append = TRUE)
  
  # Function
  cat(
    # Open chunk
    '\n```{r echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE}\n',
    '# Función para repesentar las figuras\n',
    'grid_img <- function(list_fig, dir, patt, dim = 300){
  # Abrimos todas las filas:
  forest_plots = c()
  funnel_plots = c()
  influ_plots = c()
  for (fig in list_fig){
    f1 = glue("{fig}forest.svg")
    out1 <- magick::image_read(paste0(dir, "/", f1))
    forest_plots = c(forest_plots, image_scale(out1, dim))
    
    f2 = glue("{fig}funnel.svg")
    out2 <- magick::image_read(paste0(dir,"/", f2))
    funnel_plots = c(funnel_plots, image_scale(out2, dim))  
    
    f3 = glue("{fig}influence.svg")
    out3 <- magick::image_read(paste0(dir,"/", f3))
    influ_plots = c(influ_plots, image_scale(out3, dim))
  } 
  
  lim = length(forest_plots)
  grid = NULL
  for ( i in seq(list_fig)){
    a = forest_plots[[i]]
    b = funnel_plots[[i]]
    c = influ_plots[[i]]
    
    fila <- image_append(c(a, b, c))
    
    if(is.null(grid)){
      grid = fila
    }else{
      grid = image_append(c(grid,fila),  stack = TRUE)
    } 
  }
    return(grid)
}',
    # Close chunk
    '\n```  \n\n','\n&nbsp;  \n\n',
    sep ="",
    file = report_out,
    append = TRUE)


cat(
  # Open chunk
  '\n```{r echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE, echo=FALSE, out.width=".99\\\\linewidth"}\n',
  'path = "',PlotsDir, '"\n',
  'files <- list.files(path = path, pattern = "',glue("{args$prefix}.*forest"),'")\n',
  'figuras = sapply( files, function(x) substr(x, start = 0, stop = nchar(x)-10), 
                  USE.NAMES = FALSE)\n',
  '# Plot\n',
  'grid_img(figuras, dir = path)',
  # Close chunk
  '\n```  \n\n','\n&nbsp;  \n\n',
  sep ="",
  file = report_out,
  append = TRUE)


# Create HTML
rmarkdown::render(report_out)
}


cat("\nDone\n")
cat("\n---------------------------------------------\n")

# ~~~~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~ ** ~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~~~~ #
