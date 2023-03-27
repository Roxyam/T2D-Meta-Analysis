#!/usr/bin/env Rscript

# ~~~~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~ ** ~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~~~~ #
#
# Functions.R
#
#
# Collection of udeid in this Metha-analysis
#
#
#     ****************************************************************
#     *    Author:  Roxana Andreea Moldovan                          *
#     *    Contact: roxana.andreea.moldovan@gmail.com                *
#     *    Version: 1.1                                              *
#     *    Creation date: 26/03/2023                                 *
#     ****************************************************************
#
#
# ~~~~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~ ** ~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~~~~ #


# 01Descarga

RX_GetDataGEO <- function(StudyAcc, dir=getwd()){
  "library(GEOquery)"
  setwd(dir) 
  print(glue("Processing {StudyAcc}:"))
  # Check if the file exists
  if(
    class(
      try(
        {getGEOSuppFiles(StudyAcc, makeDirectory =TRUE,
                         baseDir = dir, fetch_files = TRUE)
          # Set working directory
          setwd(glue("{dir}/{StudyAcc}"))}
      )) == "try-error"){
    # Error
    print(glue("{StudyAcc} could not be downloaded."))
  }else{
    # Decompress files
    files = list.files()
    for (file in files){
      if (endsWith(file, "tar")){
        system(glue("tar -xvf {file}"))
      } else{
        if (endsWith(file, "tar.gz")){
          system(glue("tar -xzvf {file}"))
        } else {
          if (endsWith(file, "gz")){ 
            system(glue("gzip -d {file}"))
          }
        }
      }
    }
    # Delete the compressed files
    system(paste("rm", files, sep = ' '))
    
    # Download study information
    dir.create("01RawData")
    metadata = getGEO(GEO = StudyAcc)[[1]]
    save(metadata, file = glue("{StudyAcc}_metadata.RData"))
    system( "mv `ls | grep -v 01RawData` 01RawData/")
    
    # Completed
    print(glue("Study {StudyAcc} processed successfully.."))
  }
}

RX_GetDataArrEx <- function(StudyAcc, dir=getwd()){
  # library(ArrayExpress)
  setwd(dir)
  # Creating the directory to store the data.
  path0 = glue("{StudyAcc}/01RawData")
  system(glue("mkdir {StudyAcc} {path0}"))
  
  # Data download
  if(
    class(
      try(
        {exp_set = ArrayExpress(StudyAcc, path = path0)} 
      )) == "try-error"){
    # Error
    print(glue("{StudyAcc} could not be downloaded."))
  }else{
    # Save the ExpressionSet
    save(exp_set, file = glue("{path0}/{StudyAcc}_raw.rda"))
  }
}