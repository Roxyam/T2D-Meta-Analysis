#!/usr/bin/env Rscript

# ~~~~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~ ** ~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~~~~ #
#
# GenerateAnnotPackage.R   
# 
# Generate annotation package:
#
#     This script creates a sqlite database, and then makes an
#     annotation package with it.
#
#
#     ****************************************************************
#     *    Author:  Roxana Andreea Moldovan                          *
#     *    Contact: roxana.andreea.moldovan@gmail.com                *
#     *    Version: 1.1                                              *
#     *    Creation date: 23/03/2023                                 *
#     ****************************************************************
#
#
# ~~~~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~ ** ~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~~~~ #


# ~~~~~~~~~~~~ Loading packages ~~~~~~~~~~~~ #

library(pacman)
pacman::p_load(argparse)
pacman::p_load(AnnotationForge)
pacman::p_load(glue)


# ~~~~~~~~~~~~ Parameters ~~~~~~~~~~~~ #

parser <- ArgumentParser(description="This script creates a sqlite database,
                                      and then makes an annotation package
                                      with it")

# Shemas 
schemas = available.dbschemas()
parser$add_argument("-s", "--schema",
                    action="store",
                    type="character",
                    default=NULL,
                    choices=schemas, 
                    help="Schema that you want to use to make the DB.")
# Prefix
parser$add_argument("-p", "--prefix",
                    action="store",
                    type="character",
                    default="package",
                    help="Prefix of the eventual desired package name 
                          (ie. 'prefix.db').")
# DB file
parser$add_argument("-f", "--file",
                    action="store",
                    type="character",
                    default="package",
                    required=TRUE,
                    help="The path and filename for the file to be parsed
                          and use to create the data base.")
# MapType
parser$add_argument("-m", "--maptype",
                    action="store",
                    type="character",
                    default="eg",
                    help="The type of ID that is used for the initial base
                          mapping. If using a classic base mapping file,
                          this should be the ID type present in the fileName.
                          This can be any of the following values:
                          'gb' = for genbank IDs, 
                          'ug' = unigene IDs,  
                          'eg' = Entrez Gene IDs, 
                          'refseq' = refseq IDs, 
                          'gbNRef' = mixture of genbank and refseq IDs. 
                          'eg' by default")
# Output directory
parser$add_argument("-o", "--outdir",
                    action="store",
                    type="character",
                    default=".",
                    help="Where you would like the output files to be placed,
                          by default the current directory will be taken.")
# Author
parser$add_argument("-a", "--author",
                    action="store",
                    type="character",
                    default="-",
                    help="Author or authors involved in making the package") 
# Manufacturer
parser$add_argument("-mf", "--manufacturer",
                    action="store",
                    type="character",
                    default="-",
                    help="Who made the chip being described")
# Manufacturer URL
parser$add_argument("-mu", "--manufacturerURL",
                    action="store",
                    type="character",
                    default="-",
                    help="Manufacturer website")
# Chip Name
parser$add_argument("-n", "--name",
                    action="store",
                    type="character",
                    default="-",
                    help="Name of the chip")


# ~~~~~~~~~~~~ Main ~~~~~~~~~~~~ #

#------------- Checking arguments
args <- parser$parse_args()

if (is.null(args$s)){
  cat("Schema not provided, please use uno of these:\n")
  cat(schemas, sep = ", ")
  cat("\n")
  stop("",call. = FALSE)
} else{
  cat("\nCreating SQLite based annotation package...\n\n")
}

#------------- Creation
cat("\n---------------------------------------------\n")
makeDBPackage(
  schema= args$schema,
  affy = TRUE,
  prefix = args$prefix,
  fileName = args$file, 
  baseMapType = args$maptype,
  outputDir = args$outdir,
  author = args$author,
  version = "0.1",
  manufacturer = args$manufacturer,
  manufacturerUrl = args$manufacturerURL,
  chipName = args$name
)
cat("\n---------------------------------------------\n")

# ~~~~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~ ** ~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~~~~ #
