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


# ~~~~~~~~~~~~ Functions ~~~~~~~~~~~~ #

source("Functions/Functions.R")


# ~~~~~~~~~~~~ Parameters ~~~~~~~~~~~~ #

parser <- ArgumentParser(description="This script Compute differential
                                      expression analysis of a set of
                                      studies")

# GEO studies
parser$add_argument("-g", "--geo",
                    action="store",
                    type="character",
                    default=NULL,
                    help="GEO study or studies to download
                          separated by comma.", 
                    metavar= "GEO accesion")

# Output directory
parser$add_argument("-o", "--outdir",
                    action="store",
                    type="character",
                    default=".",
                    help="Where you would like the output files to be placed,
                          by default the current directory will be taken.")


# ~~~~~~~~~~~~ Main ~~~~~~~~~~~~ #

#------------- Checking arguments
args <- parser$parse_args()

if(length(dir_in) == 1 | length(dir_in) == length(Studies)){
  files = glue("{dir_in}/{Studies}")
} else{
  stop("The input directory must be one or one for each study.")
}