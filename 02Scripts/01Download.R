#!/usr/bin/env Rscript

# ~~~~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~ ** ~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~~~~ #
#
# 01Download.R
#
#
# Data download and preparation:
#
#     This script allows you to download study data from GEO and
#     ArrayExpress using their accession numbers. It takes as input
#     the Accession numbers for GEO and/or ArrayExpress, which can be
#     entered as one or more values separated by commas. Additionally,
#     you can also specify a working directory, within which a
#     directory named Data will be created to store and process the
#     data. If no directory is specified, the current directory will
#     be used.
#
#
#     ****************************************************************
#     *    Author:  Roxana Andreea Moldovan                          *
#     *    Contact: roxana.andreea.moldovan@gmail.com                *
#     *    Version: 1.2                                              *
#     *    Creation date: 12/02/2023                                 *
#     ****************************************************************
#
#
# ~~~~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~ ** ~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~~~~ #


# ~~~~~~~~~~~~ Loading packages ~~~~~~~~~~~~ #

library(pacman)
pacman::p_load(argparse)
pacman::p_load(ArrayExpress)
pacman::p_load(GEOquery)
pacman::p_load(glue)
pacman::p_load(stringr)


# ~~~~~~~~~~~~ Functions ~~~~~~~~~~~~ #

source("Functions/Functions.R")


# ~~~~~~~~~~~~ Parameters ~~~~~~~~~~~~ #

parser <- ArgumentParser(description="Download data from GEO and
                                      BioStudies public repositories.")

# GEO studies
parser$add_argument("-g", "--geo",
                    action="store",
                    type="character",
                    default=NULL,
                    help="GEO study or studies to download
                          separated by comma.", 
                    metavar= "GEO accesion")

# ArrayExpress studies
parser$add_argument("-b", "--biostudies",
                    action="store",
                    type="character",
                    default=NULL,
                    help="BioStudies/ArrayExpress study or studies to
                          download separated by comma.", 
                    metavar= "Accesion") 

# Directory
parser$add_argument("-d", "--dir",
                    type="character",
                    default=getwd(),
                    help="Set working directory, by default the current
                          directory will be taken.",
                    metavar="path") 


# ~~~~~~~~~~~~ Main ~~~~~~~~~~~~ #

#------------- Execution
#args <- parser$parse_args(args = c('-g=c("GSE2508","GSE20950", \\
#                                   "GSE29718", "GSE64567","GSE92405",  \\
#                                   "GSE141432","GSE205668")', 
#                                   '-b="E-MEXP-1425"'))
#
#------------- Checking arguments 

if (is.null(args$g) & is.null(args$b)){
  stop("Study not specified, please enter the accession of
       one or several studies.",
       call. = FALSE)
} else{
  args$geo = eval(parse(text=(args$geo)))
  args$biostudies = eval(parse(text=(args$biostudies)))
  cat("\nDownloading data...\n\n")
}

#------------- Set and prepare directory
DirBase = args$d
DirData = glue("{DirBase}/Data")
system(glue("mkdir {DirData}"))
setwd(DirBase)

#------------- GEO
if(!is.null(args$g)){
  cat("\n\t> GEO\n")
  cat("\n---------------------------------------------\n")
  # GEO studies to download
  GEOstudies = args$g
  # Use RX function to download all stidies
  out = lapply(GEOstudies, function(e) RX_GetDataGEO(e, DirData))
  cat("\n---------------------------------------------\n")
}

#------------- BioStudies/ArrayExpress
if(!is.null(args$b)){
  cat("\n\t> BioStudies/ArrayExpress\n")
  cat("\n---------------------------------------------\n")
  # BioStudies/ArrayExpress studies to download
  Arrstudies = args$b
  # Use RX function to download all stidies
  out = lapply(Arrstudies, function(e) RX_GetDataArrEx(e, DirData))
  cat("\n---------------------------------------------\n")
}

cat("\nDone\n")
cat("\n---------------------------------------------\n")

# ~~~~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~ ** ~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~~~~ #
