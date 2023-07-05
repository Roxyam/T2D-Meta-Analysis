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
pacman::p_load(ggplot2)
pacman::p_load(ggrepel)
pacman::p_load(limma)
pacman::p_load(parallel)
pacman::p_load(stringr)
pacman::p_load(SummarizedExperiment)


# ~~~~~~~~~~~~ Functions ~~~~~~~~~~~~ #

source("Functions/Functions.R") 

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
                    help="Variable or variables to use in the
                          differential expression.")
# Covariables
parser$add_argument("-c", "--covars",
                    action="store",
                    type="character",
                    default=NULL,
                    help="Coariable or covariables to use in the
                          differential expression.")

# Data directory
parser$add_argument("-i", "--indir",
                    action="store",
                    type="character",
                    default=".",
                    help="Data directory, by default the current
                          directory will be taken.")

# Output directory
parser$add_argument("-o", "--outdir",
                    action="store",
                    type="character",
                    default=".",
                    help="Where you would like the output files to be placed,
                          by default the current directory will be taken.")

# Report
parser$add_argument("-r", "--report",
                    action="store_true",
                    help="Create an R Markdown and HTML with the results.")

# Plots
parser$add_argument("-p", "--plot",
                    action="store_true",
                    help="Add plots to report") 


# ~~~~~~~~~~~~ Main ~~~~~~~~~~~~ #

#------------- Execution
parser$print_help()
args <- parser$parse_args(args = c('-s=c("E_MEXP_1425","GSE2508","GSE20950", \\
                                   "GSE29718", "GSE64567","GSE92405",  \\
                                   "GSE141432","GSE205668")', 
                                   '-v=c("Obesity", "Diabetes")', 
                                   '-i=/home/rmoldovan/T2D-Meta-Analysis/Data',
                                   '-o=/home/rmoldovan/T2D-Meta-Analysis/04DE',
                                   '-r',
                                   '-p')) 

#------------- Prepare arguments
args$studies = eval(parse(text=(args$studies)))
args$vars = eval(parse(text=(args$vars)))
report_out = glue("{args$outdir}/DifferentialExpressionReport.rmd")

# Output directory
if(! file.exists(args$outdir)){
  system(glue("mkdir {args$outdir}"))
}

if(args$plot){
  PlotsDir = glue("{args$outdir}/Plots")
  if(! file.exists(PlotsDir)){
    # Plots directory
    system(glue("mkdir {PlotsDir}"))
  }
}


#------------- Differential expression
cat("\n\t> Differential expression\n")
cat("\n---------------------------------------------\n")
for (var in args$vars){
  Results = c()
  Studies_out = c()
  for (study in args$studies){
    Acc = study
    cat(Acc, "\n")
    DirRData = glue("{args$indir}/{Acc}/02RData")
    DirDEData = glue("{args$outdir}")
    Data = get(load(glue("{DirRData}/finalData.RData"))[[1]])
    try({      
      # Compute the differential expression
      Res <- RX_DiffExpFinal(Data, var1 = var, var2="Gender", covar = args$covars)
      # Add the annotation
      if(class(Data) == "ExpressionSet"){
        anotData = fData(Data)  
      }else{
        anotData = rowData(Data)      
      }
      
      Res = sapply(names(Res), 
                   function(i) setNames(mclapply(names(Res[[i]]),
                                        function(j){ # Top table
                                          Res[[i]][[j]]$"TopTab" = cbind(Res[[i]][[j]]$"TopTab",
                                                         anotData[match(rownames(Res[[i]][[j]]$"TopTab"),
                                                                        anotData$ENTREZID),])
                                          # Volcano plot
                                          if(isTRUE(args$plot)){
                                            Plot = RX_VolcanoPlot(Res[[i]][[j]]$"TopTab", logFC_lim = 0) 
                                            ggsave(Plot,
                                                   device = "svg",
                                                   height = 6, 
                                                   width = 10,
                                                   filename=str_replace_all(glue("{DirDEData}/Plots/{study}_{i}_{j}.svg"), "\\s", ""))
                                          }
                                        #Res[[i]][[j]]$"Plot" = RX_VolcanoPlot(Res[[i]][[j]]$"TopTab", logFC_lim = 0)
                                        return(Res[[i]][[j]])} #, simplify = FALSE
                                      ), nm=names(Res[[i]])) ,simplify = FALSE
                   )
      
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
pacman::p_load(ggpubr)
pacman::p_load(flextable)
pacman::p_load(officer)
pacman::p_load(dplyr)
pacman::p_load(magick)

# Functions
#source("Functions/Functions.R")
Dir = "',args$outdir,'" ',
    # Close chunk
    '\n```  \n\n','\n&nbsp;  \n\n',
    sep ="",
    file = report_out,
    append = TRUE) 
  # Load packages
  cat(
    # Open chunk
    '\n```{r echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE}\n',
    '# Plot output
RX_grid_img <- function(paths){
  grid = NULL
  for (i in seq(from = 1, to=length(paths), by=2)){
    a = image = image_read(paths[i])
    b = image_read(paths[i+1])
    row = image_append(c(a, b))
    if(is.null(grid)){
      grid = row
    }else{
      grid = image_append(c(grid,row),  stack = TRUE)
    } }
    return(grid)
    }\n',
    '# Output table
RX_table <- function(Data,
                     Color1 = "#007FA6",
                     Color2 = "#007FA6" ,
                     Color.l = "black" ,
                     Color.l2 = "#878787",
                     Color.up = "#CB326D",
                     Color.down = "#00AFBB",
                     footer = NULL,
                     font = "Calibri"){
    
  cols=colnames(Data)[-c(1,2)] # Columns with info
  col_keys = colnames(Data) # Nombres unicos de las columnas
  # Header
  header = data.frame(col_keys = col_keys,
                      H1 = c("Contraste", "Contraste", rep("Estudios", length(cols))),
                      H2 = c("Contraste", "Contraste", cols)
                      )
  
  ## Flextable
  # Set info
  ft <- flextable(Data, col_keys = colnames(Data))
  # Set header
  ft <- set_header_df(ft, mapping = header, key = "col_keys" )
  # Format number
  ft <- colformat_num(ft, big.mark = ".", decimal.mark = "," )
  # Set footer
  if(! is.null(footer)){
    ft <- add_footer(ft, values = footer) %>% 
      merge_h(part = "footer")  %>% 
      bold(j=1, part = "footer") %>%
      hline(border = fp_border(width = 2, color = Color1), part = "footer")
  }
  
  # Merge headers
  ft <- merge_h(ft, part = "header") %>% 
    merge_v( part = "header")
  
  # Merge contrastes
  ft <- merge_v(ft, j = c(1,2), part = "body")
  
  # Horizontal lines
  ft <- hline(ft, i=seq(from=3, to=nrow(Data), by=3),
              border = fp_border(width = 1, color = Color.l2), part="body")
  ft <- hline(ft,border = fp_border(width = 2, color = Color1), part = "header")
  
  # Outer bottom
  ft <- hline_bottom(ft, border = fp_border(width = 2, color = Color1))
  
  # Alineamiento
  ft <- flextable::align(ft, align = "center", part = "header")
  ft <- flextable::align(ft, align = "right", part = "footer")
  ft <- flextable::align(ft, j=1, align = "left", part = "footer")
  
  
  # Header design
  ft <- fontsize(ft, i = 1:2, size = 12, part = "header") %>%
    color(i=1, color = Color2, part = "header") %>%
    color(i=2, color = Color.l2, part = "header") %>%
    bold(i=1, bold = TRUE, part = "header")
  # Rotate names
  #ft <- rotate(ft, i = 2, rotation = "btlr", part = "header", align = "bottom")
  
  # Font
  #font(ft, fontname = font, part = "all")
  
  # Matriz de colores
  colormatrix <- ifelse(Data[, cols] == 0, Color.l2, Color.l )
  # Color
  ft <- color(ft, j=cols, color = colormatrix, part = "body") 
  # Color
  x = Data[, 2]
  colors = case_when(
    x == "Up" ~ Color.up,
    x == "Down" ~ Color.down,
    TRUE ~ Color.l
  )
  ft <- color(ft, j=2, color = colors, part = "body") %>%
    bold(j=2, bold = TRUE, part = "body")
  return(ft)
    
}',
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
      'Res = get(load(glue("{Dir}/DifferentialExpression', var,'.RData")))\n',
      'St = names(Res)\n',
      'Toptabs = sapply(names(Res), 
                 function(est) sapply(names(Res[[est]]), 
                                      function(tis) sapply(Res[[est]][[tis]], 
                                                           function(x) x$TopTab, 
                                                           simplify = FALSE),
                                      simplify = FALSE),
                 simplify = FALSE) ', 
      '\nPlots = sapply(names(Res), 
                 function(est) sapply(names(Res[[est]]), 
                                      function(tis) sapply(Res[[est]][[tis]], 
                                                           function(x) x$Plot, 
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
      # Number of tissues
      nSt = sum(sapply(Results, function(x)  ! is.null(x[[tis]])))
      cat(
        # Title
        '\n#### ',tis,'  {.tabset .tabset-fade -}  \n',
        '\n&nbsp;  \n',
        '\n##### Tabla  {-}  \n\n',
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
                              function(x1){ up = length(which(x1$logFC > 0))
                                            down = length(which(x1$logFC < 0))
                                            return(c(up,down, sum(up,down)))},
                              simplify = FALSE, USE.NAMES = TRUE),
           simplify = FALSE, USE.NAMES = TRUE)

rows = unlist(sapply(Resume_tis, function(x) names(x), USE.NAMES = TRUE, simplify = FALSE))
rows = unique(rows)
# Sep
cols = names(Resume_tis)

DT_tis0 = sapply(cols,
                 function(col) sapply(rows,
                                      function(row) RX_ifelse(is.null(Resume_tis[[col]][[row]]), rep("-", 3), Resume_tis[[col]][[row]]),
                                      USE.NAMES = T),
                 USE.NAMES = T)
DT_tis = data.frame("Contrast" = rep(rows, each = 3),
                    "Direction" = rep(c("Up", "Down", "Total"), times = length(rows)),
                    DT_tis0)
#datatable(DT_tis0, caption = "Resultados expresión diferencial ', tis,'")
foot = sapply(names(Toptabs_tis), function(x) nrow(Toptabs_tis[[x]][[1]]))
val = c("Nº total de genes", "Nº total de genes",foot)
names(val)[1:2] = c("Contrast", "Direction")
RX_table(DT_tis, footer = val) # OUT', 
        # Close chunk
        '\n```  \n\n','\n&nbsp;  \n\n',
        '\n***  \n',
        sep ="",
        file = report_out,
        append = TRUE)
      
      # Volcano plot
  if(args$plot){      
  cat(
        '\n##### Plots  {.tabset .tabset-fade -}  \n\n',
        sep ="",
        file = report_out,
        append = TRUE)
      for (St in args$studies){
        
        # Recorremos los contrastes por tejido
        contrasts = names(Results[[St]][[tis]]) 
        if(!is.null(contrasts)){
        cat(
          '\n###### ', St ,'  {-}  \n\n',
          # Open chunk
          '\n```{r echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE, fig.width=12, fig.height=', (ceiling(length(contrasts)/2)*5),', fig.align="center"}\n\n',
          'contrasts = c(', paste("'",str_replace_all(contrasts,"\\s", ""),"'", collapse = ", ", sep = ""),')\n',
          'paths = glue("{Dir}/Plots/',St,'_',tis,'_{contrasts}.svg")\n',
          'RX_grid_img(paths)\n', 
          # Close chunk
          '\n```  \n\n','\n&nbsp;  \n\n',
          '\n***  \n',
          sep ="",
          file = report_out,
          append = TRUE)
          }
        
      }}
    }
  }
# Create HTML
rmarkdown::render(report_out)
}

cat("\nDone\n")
cat("\n---------------------------------------------\n")
