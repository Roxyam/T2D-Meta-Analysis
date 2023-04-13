#!/usr/bin/env Rscript

# ~~~~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~ ** ~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~~~~ #
#
# 05MetaAnalysis.R
#
#
# Computes the meta-analysis for a given set of studies:
#
#     
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

source("../Functions/Functions.R")

RX_getSE <- function(Data # Fit limma
){
  # Método 1:
  TopTab <- topTable(Data, number = "all", confint=TRUE, adjust.method = "fdr", sort.by = "none")
  # Calculamos el SE
  TopTab[, "SE"] <- (TopTab[, "CI.R"] - TopTab[, "CI.L"])/ 3.92
  
  # Método 2:
  coefSE = sqrt(Data$s2.post) * Data$stdev.unscaled
  TopTab[, c("coef", "coefSE")] = c(Data$coefficients, coefSE)
  
  return(TopTab)
}

RX_MetaAnalysis_prep <- function(Studies, # Studies accesion to include
                                 Tissue = c("SAT", "VAT", "LSAT"), # One or several tissues to include in the meta-analysis to include
                                 Contrast, # Contrast to use
                                 Datas, # data
                                 SE_coef = c("SE", "coefSE") # SE to take
){
  # Check if the arguments are suitable
  match.arg(Tissue, several.ok = TRUE)
  #match.arg(Disease, several.ok = FALSE)
  match.arg(SE_coef, several.ok = FALSE)
  
  # Extract the information from the indicated tissues and contrast
  Data0 = sapply(Studies,
                 function(st) sapply(Tissue,
                                     function(ts) sapply(Contrast,
                                                         function(c) Datas[[st]][[ts]][[c]], 
                                                         simplify = FALSE),
                                     simplify = FALSE),
                 simplify = FALSE)
  
  Data = unlist(unlist(Data0, recursive = FALSE), recursive = FALSE)
  names(Data) = Studies
  Data = Filter(Negate(is.null), Data)
  
  #------------- 1. Preparing input for meta-analysis
  # Get SE
  DataSE = sapply(names(Data), function(x) RX_getSE(Data[[x]]), simplify = FALSE) 
  # Get Gene and annotation
  IDs = unlist(sapply(DataSE, function(x) rownames(x), simplify = TRUE), use.names = FALSE)
  SYMBOL = sapply(Data, function(x) x$TopTab[,c("ENTREZID", "ENSEMBL", "SYMBOL")], simplify = FALSE)
  ANOT = do.call("rbind", SYMBOL)
  cat(glue("\n{length(unique(IDs))} of the {length(IDs)} genes have been taken."), "\n") 
  IDs = unique(IDs)
  ANOT = ANOT[match(unique(IDs), ANOT$ENTREZID),]
  rownames(ANOT) = ANOT$ENTREZID
  
  # LogFC matrix
  MatLogFC <- sapply(DataSE, function(x) x[IDs, "logFC"]) 
  rownames(MatLogFC) <- IDs
  
  # SD matrix
  MatSE <- sapply(DataSE, function(x) x[IDs, SE_coef]) 
  rownames(MatSE) <- IDs
  
  # Filter keeping the genes that appear in more than one study
  keep = which(rowSums(!is.na(MatSE)) > 1)
  MatLogFC = MatLogFC[keep,]
  MatSE = MatSE[keep,] 
  
  return(list("MatLogFC" = MatLogFC, "MatSE" = MatSE, Anot = ANOT))
}



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
                    help="Contrast or contrasts to include.")

# Out dir
parser$add_argument("-o", "--outdir",
                    action="store",
                    type="character",
                    default=".",
                    help="Where you would like the output files to be placed,
                          by default the current directory will be taken.")

# Output files prefix
parser$add_argument("-p", "--prefix",
                    action="store",
                    type="character",
                    default="Meta-analysis",
                    required=TRUE,
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
                    dafault=TRUE,
                    help="Create an R Markdown and HTML with the results.")
# Repeat
parser$add_argument("-re", "--redo",
                    action="store_true",
                    default = FALSE,
                    help="Repeat the analysis, if it is not specified,
                          the files of previous executions are taken")

# ~~~~~~~~~~~~ Main ~~~~~~~~~~~~ #

#------------- Checking arguments
args <- parser$parse_args()

### XX NONO
args = list()
args$studies = "ALL"
args$tissue = "SAT"
args$filein = "../Data/DiffExprs_Obesity.RData"
args$outdir = "."
args$prefix = "Meta-analysis"
args$contrast = "Ob - C, Ob.M - C.M, Ob.F - C.F, (Ob.M - C.M) - (Ob.F - C.F)"
args$report = TRUE
args$redo = FALSE

#load(file="../Data/DE/DifferentialExpressionObesity.RData")
#Rscript 05MetaAnalysis.R -t SAT -f "../Data/DiffExprsObesity.RData" -o ../Data/Meta-analysis_Obesity_all.RData -c "Ob - C, Ob.M - C.M, Ob.F - C.F, (Ob.M - C.M) - (Ob.F - C.F)"


# Load the data
Datas = get(load(args$filein))

# Get the studies
if (args$studies != "ALL"){
  # Take studies to include
  Studies = stringr::str_split_1(args$studies, pattern = ",")
  Studies = trimws(Studies, whitespace = "[ \t\r\n\\.]")
} else{
  Studies = names(Datas)
}

# Cheek
if (! all(Studies %in% names(Datas))){
  cat("The indicated file does not contain information on",
      "the following studies. \n")
  cat("\t", Studies[! Studies %in% names(Datas)])
  cat("\n")
  stop("Please review input.",call. = FALSE)
}

# Get the tissues
Tissue = stringr::str_split_1(args$tissue, pattern = ",")
Tissue = trimws(Tissue, whitespace = "[ \t\r\n\\.]")

# Get the contrasts
Contrast = stringr::str_split_1(args$contrast, pattern = ",")
Contrast = trimws(Contrast, whitespace = "[ \t\r\n\\.]")

#------------- 1. Meta-analysis for genes
cat("\n---------------------------------------------\n")
cat("Analyzing...\n")

resultMAs = mclapply(Contrast, 
                     function(C) RX_MetaAnalysis_prep(Studies = Studies,
                                                      Tissue = Tissue,
                                                      Datas = Datas, 
                                                      Contrast = C,
                                                      SE_coef = args$method))
names(resultMAs) = Contrast

# Save
resultMAs = list()
save(resultMAs, file = glue("{args$outdir}/{args$prefix}.RData")) 

cat("\n---------------------------------------------\n")

