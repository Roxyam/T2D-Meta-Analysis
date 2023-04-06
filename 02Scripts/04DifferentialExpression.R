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

#------------- Checking arguments
args <- parser$parse_args()

