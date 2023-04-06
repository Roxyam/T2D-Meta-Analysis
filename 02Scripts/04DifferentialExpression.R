#!/usr/bin/env Rscript

# ~~~~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~ ** ~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~~~~ #
#
# 04DifferentialExpression.R
#
#
# Compute differential expression analysis of a set of studies:
#
#     This script takes the data referring to a set of studies and
#     computes the differential expression analysis independently 
#     for each one of them, carrying out the contrasts for one or a
#     set of variables and its intersection with the sex of the
#     patient. 
#     In addition to saving the results obtained in RData format,
#     create a report that shows them.
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
pacman::p_load(Biobase)
pacman::p_load(DT)
pacman::p_load(glue)
pacman::p_load(limma)
pacman::p_load(stringr)
pacman::p_load(SummarizedExperiment) 


# ~~~~~~~~~~~~ Functions ~~~~~~~~~~~~ #

source("Functions/Functions.R")
RX_contrast <- function(expmatrix, design, C){
  # library(limma)
  set.seed(2808) 
  contMatrix <- makeContrasts(contrasts = C, levels = design)
  fit <- lmFit(expmatrix, design)
  fit2 <- contrasts.fit(fit, contMatrix)
  Res <- eBayes(fit2)
  Res$TopTab <- topTable(Res, number = Inf, adjust.method = "BH") #, sort.by = "logFC"
  
  return(Res)
}

RX_DiffExp <- function(phenoData, expressData , C, var="Group", studyType = "Array"){
  #library(Biobase)
  #library(stringr)
  set.seed(2808)
  if(length(var) == 2){
    contrast = interaction(phenoData[,var[1]], phenoData[,var[2]])
  }else{
    contrast = phenoData[,var]
  }
  
  design = model.matrix(~0+contrast)
  colnames(design) <- str_replace(levels(contrast) , "-", "_")
  
  # Si es RNA-seq usamos voom
  if (studyType == "RNA-seq"){
    expressData = voom(expressData, design = design, plot = F)
  }
  # Top table
  Res <- sapply(C, function(x) RX_contrast(expressData, design, x),
                simplify = FALSE, USE.NAMES = TRUE)
  return(Res)
}

# Global
RX_DiffExpFinal <- function(Data, var1="Group", var2="Gender"){
  # library(SummarizedExperiment)
  # library(Biobase)
  # library(glue) 
  set.seed(2808) 
  if(class(Data) == "ExpressionSet"){
    phenoData = pData(Data)
    expressData = exprs(Data)
    studyType = "Array"
  }else{
    phenoData = colData(Data)
    studyType = "RNA-seq"
    expressData = assay(Data)      
  }
  # Tomamos las variables
  Group = str_replace(levels(phenoData[,var1]) , "-", "_")
  Gender = levels(phenoData[,var2])
  Tissue = levels(phenoData$Tissue)
  # Separamos por tejidos
  phenoDatas = sapply(Tissue, function(x) phenoData[which(phenoData$Tissue == x),],
                      simplify = FALSE, USE.NAMES = TRUE)
  expressDatas = sapply(Tissue, function(x) expressData[,which(phenoData$Tissue == x)],
                        simplify = FALSE, USE.NAMES = TRUE)
  
  V1 = str_replace(levels(phenoData[,var1]) , "-", "_")
  V2 = levels(phenoData[,var2])
  # Contrastes
  'C1 = c(glue("{V1[2]} - {V1[1]}"))
  C2 = c(glue("({Group[2]}.{Gender} - {Group[1]}.{Gender})"))
  C3 = c(C2, glue("{C2[1]} - {C2[2]}"))'
  # Contrastes
  C1 = combn(rev(Group), m=2, FUN = paste,collapse=" - ", simplify = TRUE)
  
  C2 = sapply(Gender, 
              function(x) combn(levels(interaction(Group, x))[match(rev(Group), strsplit2(levels(interaction(Group, x)),
                                                                                          split = ".", fixed = TRUE)[,1])],
                                m=2, FUN = paste,collapse=" - ",
                                simplify = TRUE), simplify = F)
  
  C3 = sapply(seq(C1), function(i) glue("({C2$M[i]}) - ({C2$F[i]})"))
  C3 = c(unlist(C2), C3, use.names = FALSE)
  
  
  Tt = sapply(Tissue, function(x) c(RX_DiffExp(phenoDatas[[x]], expressDatas[[x]] , C1, var=var1, studyType = studyType),
                                    RX_DiffExp(phenoDatas[[x]], expressDatas[[x]], C3, var=c(var1, var2), studyType = studyType)),
              simplify = FALSE, USE.NAMES = TRUE) 
  
  return(Tt)
}


# ~~~~~~~~~~~~ Parameters ~~~~~~~~~~~~ #

parser <- ArgumentParser(description="This script Compute differential
                                      expression analysis of a set of
                                      studies")

# Studies
parser$add_argument("-s", "--studies",
                    action="store",
                    type="character",
                    required=TRUE,
                    help="Studies to include separated by comma")

# Variables
parser$add_argument("-v", "--vars",
                    action="store",
                    type="character",
                    default="Group",
                    choices = c("Group","Obesity", "Diabetes"), 
                    help="Variable or variables to use in the
                          differential expression.")

# Data directory
parser$add_argument("-o", "--outdir",
                    action="store",
                    type="character",
                    default=".",
                    help="Where you would like the output files to be placed,
                          by default the current directory will be taken.")

# Output directory
parser$add_argument("-i", "--indir",
                    action="store",
                    type="character",
                    default=".",
                    help="Data directory, by default the current
                          directory will be taken.")

# Report
parser$add_argument("-r", "--report",
                    action="store",
                    type="character",
                    dafault=TRUE,
                    help="Create an R Markdown and HTML with the results.")


# ~~~~~~~~~~~~ Main ~~~~~~~~~~~~ #

#------------- Prepare arguments
args <- parser$parse_args()
report_out = glue("{args$outdir}/DifferentialExpressionReport.rmd")

#------------- Differential expression
cat("\n\t> Differential expression\n")
cat("\n---------------------------------------------\n")
for (var in args$vars){
  Results = c()
  Studies_out = c()
  for (study in args$studies){
    Acc = study
    DirRData = glue("{args$indir}/{Acc}/02RData")
    DirDEData = glue("{args$outdir}")
    load(glue("{DirRData}/finalData.RData"))
    Data = finalData
    try({      
      # Compute the differential expression
      Res <- RX_DiffExpFinal(Data, var1 = var, var2="Gender")
      # Add the annotation
      if(class(Data) == "ExpressionSet"){
        anotData = fData(Data)  
      }else{
        anotData = rowData(Data)      
      }
      
      Res = sapply(names(Res), 
                   function(i) sapply(names(Res[[i]]),
                                      function(j){ cbind(Res[[i]][[j]]$"TopTab",
                                                         anotData[match(rownames(Res[[i]][[j]]$"TopTab"),
                                                                        anotData$ENTREZID),])
                                        return(Res[[i]][[j]])},
                                      simplify = FALSE),
                   simplify = FALSE)  
      # Save the information in a R data file
      Results = c(Results, list(Res))
      Studies_out = c(Studies_out, study)
      #save(Res, file = glue("{DirDEData}/{Acc}DifferentialExpression{var}.RData")) 
    }, silent = TRUE)
  }
  names(Results) = Studies_out
  save(Results, file = glue("{DirDEData}/DifferentialExpression{var}.RData")) 
}

#------------- Report
cat("\n\t> Report\n")
cat("\n---------------------------------------------\n")
if(isTRUE(args$report)){
  # File info
  cat('---
title: "04 Differential expression analysis"
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
  
  ## Análisis de expresión diferencial {.tabset .tabset-fade -} 
  cat('\n## Análisis de expresión diferencial {.tabset .tabset-fade -}  \n',
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
# Functions
source("../Functions/Functions.R")',
    # Close chunk
    '\n```  \n\n','\n&nbsp;  \n\n',
    sep ="",
    file = report_out,
    append = TRUE)
  
  # Information by comparison
  for (var in args$vars){
    Results = get(load(glue("{args$outdir}/DifferentialExpression{var}.RData")))
    cat('\n### ', var ,'  \n',
        '&nbsp;  \n\n',
        sep = "",
        file = report_out,
        append = TRUE)
    # Read data
    cat(
      # Open chunk
      '\n```{r echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE}\n',
      '# Read data\n',
      'Res = get(load("',args$outdir,'/DifferentialExpression', var,'.RData"))\n',
      'St = names(Res)\n',
      'Toptabs = sapply(names(Res), 
                 function(est) sapply(names(Res[[est]]), 
                                      function(tis) sapply(Res[[est]][[tis]], 
                                                           function(x) x$TopTab, 
                                                           simplify = FALSE),
                                      simplify = FALSE),
                 simplify = FALSE) ',  
      # Close chunk
      '\n```  \n\n','\n&nbsp;  \n\n',
      sep ="",
      file = report_out,
      append = TRUE)
    
    # Take all tissues
    Tissues = unique(unlist(sapply(Results, function(x) names(x)), use.names = FALSE))
    for (tis in Tissues){
      cat(
        # Title
        '\n#### ',tis,'  \n',
        '\n&nbsp;  \n',
        # Open chunk
        '\n```{r echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE}\n\n',
        'Toptabs_tis = sapply(Toptabs, 
           function(x) x[names(x) == "', tis,'"],
           simplify = FALSE, USE.NAMES = TRUE) 

Simp = sapply(Toptabs_tis, function(x) RX_ifelse(is.null(x$', tis,'[1]), NULL, x$', tis,'),
              simplify = FALSE) 
Toptabs_tis = Filter(Negate(is.null), Simp)

Toptabs_tis_sig = sapply(Toptabs_tis, 
           function(x) sapply(x, 
                              function(x1) x1[x1$"adj.P.Val" <= 0.05,],
                              simplify = FALSE, USE.NAMES = TRUE),
           simplify = FALSE, USE.NAMES = TRUE) 

Resume_tis = sapply(Toptabs_tis_sig, 
           function(x) sapply(x, 
                              function(x1) nrow(x1),
                              simplify = FALSE, USE.NAMES = TRUE),
           simplify = FALSE, USE.NAMES = TRUE)

rows = unlist(sapply(Resume_tis, function(x) names(x), USE.NAMES = TRUE, simplify = FALSE))
rows = unique(rows)
# Por separado
cols = names(Resume_tis)
DT_tis0 = sapply(cols, function(col) sapply(rows, function(row) RX_ifelse(is.null(Resume_tis[[col]][[row]]), "-", Resume_tis[[col]][[row]]), USE.NAMES = T), USE.NAMES = T)
datatable(DT_tis0, caption = "Resultados expresión diferencial ', tis,'")

# Conjunto 
cols= St
DT_tis = sapply(cols, function(col) sapply(rows, function(row) RX_ifelse(is.null(Resume_tis[[col]][[row]]), "-", Resume_tis[[col]][[row]]), USE.NAMES = T), USE.NAMES = T)
DT_tis = datatable(DT_tis, caption = "Resultados expresión diferencial ', tis,'")', 
        # Close chunk
        '\n```  \n\n','\n&nbsp;  \n\n',
        '\n***  \n',
        sep ="",
        file = report_out,
        append = TRUE)
    }
  }
}

# Create HTML
rmarkdown::render(report_out)

cat("\nDone\n")
cat("\n---------------------------------------------\n")
