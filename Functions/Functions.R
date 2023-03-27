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

RX_GetDataGEO <- function(estudio, dir=getwd()){
  "library(GEOquery)"
  setwd(dir) 
  print(glue("Processing {estudio}:"))
  # Check if the file exists
  if(
    class(
      try(
        {getGEOSuppFiles(estudio, makeDirectory =TRUE,
                         baseDir = dir, fetch_files = TRUE)
          # Set working directory
          setwd(glue("{dir}/{estudio}"))}
      )) == "try-error"){
    # Error
    print(glue("{estudio} could not be downloaded."))
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
    metadata = getGEO(GEO = estudio)[[1]]
    save(metadata, file = glue("{estudio}_metadata.RData"))
    system( "mv `ls | grep -v 01RawData` 01RawData/")
    
    # Completed
    print(glue("Study {estudio} processed successfully.."))
  }
}

RX_GetDataArrEx <- function(estudio, dir=getwd()){
  # library(ArrayExpress)
  setwd(dir)
  # Creamos el directorio para almacenar los datos
  path0 = glue("{estudio}/01RawData")
  system(glue("mkdir {estudio} {path0}"))
  #system(paste0("cd ", estudio))
  # Descargamos los datos
  if(
    class(
      try(
        {exp_set = ArrayExpress(estudio, path = path0)} 
      )) == "try-error"){
    # Error
    print(glue("{estudio} no ha podido ser descargado."))
  }else{
    # Guardado del ExpressionSet
    save(exp_set, file = glue("{path0}/{estudio}_raw.rda"))
  }
}