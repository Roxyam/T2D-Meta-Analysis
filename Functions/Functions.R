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


# ~~~~~~~~~~~~ 01Download ~~~~~~~~~~~~ #

RX_GetDataGEO <- function(StudyAcc, # GEO Accession
                          dir=getwd() # Working directory, by default the current directory is taken
                          ){
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

RX_GetDataArrEx <- function(StudyAcc, # ArrayExpress Accession
                            dir=getwd() # Working directory, by default the current directory is taken
                            ){
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


# ~~~~~~~~~~~~ 02Construct ~~~~~~~~~~~~ #

RX_probe2gene <- function(Data, # ExpressionSet or SummarizedExperiment from microarray study
                          anotdb = org.Hs.eg.db, # Pre-Loaded Annotation Package
                          from = "PROBEID", # Current identifiers
                          to = "ENTREZID", # Identifiers to take
                          multi.to = "first", # Treatment of multiple matches with destination identifiers
                          multi.from = median # Treatment of multiple matches with source identifiers. It corresponds to a specific function, by default the median is taken.
                          ){  
  # library(AnnotationDbi)
  if (class(Data) == "SummarizedExperiment"){ 
    DataOut = Data
    exData = assay(Data)
  }else{
    DataOut = Data 
    exData = exprs(Data)
    #rowData(DataOut) = data.frame(rowData(DataOut))
  }
  # Annotation
  Anot0 = AnnotationDbi::select(x = anotdb,
                                keys = rownames(DataOut),
                                columns = to,
                                keytype = from)
  
  
  # If there is more than one gene per probe, we take the first one.
  ind1 = match(rownames(DataOut), Anot0[,from])
  Anot = Anot0[ind1,]
  
  # Delete probes that do not correspond to any gene. 
  Anot = Anot[which(! is.na(Anot[,to])),]
  
  # Join the annotation with the input data
  Exp0 = data.frame(from = rownames(exData), exData)
  Exp1 = merge(Anot, Exp0, by.x = glue("{from}"), by.y = 1)
  
  # Take the samples name
  Samp = colnames(DataOut)
  
  # If we have more than one probe per gene we use the 'multi.from' function.
  Exp2 = aggregate(Exp1[which(colnames(Exp1) %in% Samp)], 
                   by = list(X= Exp1[,to]),
                   multi.from)
  
  # We have in X the genes according to the indicated nomenclature.
  # Order as in the original,
  ord = match(unique(Exp1[,to]),  Exp2$X)
  Exp3 <- Exp2[ord,]
  ind2 = match(Exp3[,"X"], Anot[,to])
  Anot2 = Anot[ind2,] # Anot 2 es la anotacion final
  
  # Match the data
  ind3 = match(rownames(DataOut), Anot2[,from])
  # Assign the values to the original 
  if (class(Data) == "SummarizedExperiment"){ 
    rowData(DataOut) = Anot2[ind3,]
    DataOut <- DataOut[which(!is.na(rowData(DataOut)$PROBEID),)]
    rownames(DataOut) = Exp3$X
    rownames(Exp3) = Exp3$X
    # We assign the values to the array
    assay(DataOut, withDimnames=FALSE) <- as.matrix(Exp3[,-1])
  }else{
    fData(DataOut) = (Anot2[ind3,])
    DataOut <- DataOut[which(!is.na(fData(DataOut)$PROBEID),)]
    rownames(DataOut) = Exp3$X
    rownames(Exp3) = Exp3$X
    # We assign the values to the array
    exprs(DataOut) <- as.matrix(Exp3[,-1]) 
  }
  # Return
  return(DataOut)
}


RX_annot <- function(Data, # ExpressionSet or SummarizedExperiment 
                     db = org.Hs.eg.db,  # Pre-Loaded Annotation Package
                     fromAnot="ENTREZID",  # Current identifiers
                     toAnot = c("ENSEMBL","SYMBOL") # Identifiers to take
){
  #library("AnnotationDbi")
  #library("BiocGenerics")
  #library(org.Hs.eg.db)
  Dataout = Data
  # Search for matches
  anot = AnnotationDbi::select(x = db,
                               keys = rownames(Dataout),
                               columns = toAnot,
                               keytype = fromAnot,
                               multiVals = "first")
  
  # Check and deal with multiple matches
  if(! isTRUE(all.equal(rownames(Dataout),anot[,fromAnot]))){  
    ind = BiocGenerics::match(rownames(Dataout), anot[,fromAnot]) 
    if (class(Data) == "SummarizedExperiment"){
      # Combine with the above information
      rowData(Dataout) =cbind(anot[ind,],
                              rowData(Dataout)[which(!colnames(rowData(Dataout)) 
                                                     %in% c(fromAnot,toAnot))])
    }else{
      # Combine with the above information
      fData(Dataout) =cbind(anot[ind,],
                            fData(Dataout)[which(!colnames(fData(Dataout)) 
                                                 %in% c(fromAnot, toAnot))])
      
    }
  }
  
  return(Dataout)
}

