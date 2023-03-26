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

## Funciones empleadas
RX_GetDataGEO <- function(estudio, dir=getwd()){
  "library(GEOquery)"
  setwd(dir) 
  print(glue("Procesando {estudio}:"))
  # Comprobamos si los fichero existen
  if(
    class(
      try(
        {getGEOSuppFiles(estudio, makeDirectory =TRUE,
                         baseDir = dir, fetch_files = TRUE)
          # Descompresion
          setwd(glue("{dir}/{estudio}"))}
      )) == "try-error"){
    # Error
    print(glue("{estudio} no ha podido ser descargado."))
  }else{
    # Descomprimimos solo los raw
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
    #Eliminar los comprimidos
    system(paste("rm", files, sep = ' '))
    
    # Descargamos la información del estudio
    dir.create("01RawData")
    metadata = getGEO(GEO = estudio)[[1]]
    save(metadata, file = glue("{estudio}_metadata.RData"))
    system( "mv `ls | grep -v 01RawData` 01RawData/")
    
    # Completado
    print(glue("El estudio {estudio} se ha procesado correctamente."))
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